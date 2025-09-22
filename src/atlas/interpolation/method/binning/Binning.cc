/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/binning/Binning.h"

#include <climits>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/library/config.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"
#include "atlas/mesh/actions/GetCubedSphereNodalArea.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "eckit/config/LocalConfiguration.h"


namespace atlas {
namespace interpolation {
namespace method {

namespace {

MethodBuilder<Binning> __builder("binning");

using TripletType = linalg::Triplet<Binning::ValueType, Binning::IndexType>;
static_assert(std::is_trivially_copyable_v<TripletType>);

class MpiBuffer {
public:
    MpiBuffer(): buffer_{mpi::comm().size()} {}

    std::size_t size(std::size_t rank) const {
        const auto n = buffer_.at(rank).size();
        ATLAS_ASSERT(n % sizeof(TripletType) == 0, "Buffer size not a multiple of Triplet size");
        return n / sizeof(TripletType);
    }

    void pushBack(std::size_t rank, const TripletType& triplet) {
        auto& subBuffer                        = buffer_.at(rank);
        alignas(TripletType) auto appendBuffer = std::array<char, sizeof(TripletType)>{};
        std::memcpy(appendBuffer.data(), &triplet, sizeof(TripletType));
        subBuffer.insert(subBuffer.end(), appendBuffer.begin(), appendBuffer.end());
    }

    auto get(std::size_t rank) const {
        const auto& subBuffer = buffer_.at(rank);
        return [&subBuffer](size_t index) -> TripletType {
            TripletType triplet{};
            std::memcpy(&triplet, subBuffer.data() + index * sizeof(TripletType), sizeof(TripletType));
            return triplet;
        };
    }

    MpiBuffer allToAll() const {
        auto recvBuffer = MpiBuffer{};
        mpi::comm().allToAll(buffer_, recvBuffer.buffer_);
        return recvBuffer;
    }

private:
    std::vector<std::vector<char>> buffer_{};
};

}  // namespace

Binning::Binning(const Config& config): Method(config) {
    const auto* conf = dynamic_cast<const eckit::LocalConfiguration*>(&config);
    ATLAS_ASSERT(conf, "config must be derived from eckit::LocalConfiguration");
    interpAncillaryScheme_ = conf->getSubConfiguration("scheme");
}


void Binning::do_setup(const Grid& source, const Grid& target, const Cache&) {
    ATLAS_NOTIMPLEMENTED;
}

void Binning::do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) {
    ATLAS_NOTIMPLEMENTED;
}


void Binning::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    ATLAS_TRACE("atlas::interpolation::method::Binning::do_setup()");

    source_ = source;
    target_ = target;

    if (target_.size() == 0) {
        return;
    }

    const auto interp              = Interpolation(interpAncillaryScheme_, target_, source_);
    const auto interpMatrixStorage = interpolation::MatrixCache(interp).matrix();
    const auto interpMatrixView    = linalg::make_host_view<ValueType, IndexType>(interpMatrixStorage);
    auto transposeMatrixStorage    = transposeAndHaloExchange(interpMatrixView);
    auto transposeMatrixView       = linalg::make_host_view<ValueType, IndexType>(transposeMatrixStorage);

    const IndexType transposeRows = transposeMatrixView.rows();
    const IndexType transposeCols = transposeMatrixView.cols();


    const auto areaWeightsVector = getAreaWeights();

    auto triplets = std::vector<TripletType>{};
    triplets.reserve(transposeMatrixView.nnz());

    const auto weightedValue = [&](IndexType col, ValueType value) {
            return value * areaWeightsVector.at(col);
    };

    // Area-weight and normalise rows.
    for (std::size_t rowIdx = 0; rowIdx < transposeRows; ++rowIdx) {

        ValueType rowSum = 0.;
        linalg::sparse_matrix_for_each_triplet(rowIdx, transposeMatrixView,
                                               [&](IndexType row, IndexType col, ValueType value) { rowSum += weightedValue(col, value); });

        linalg::sparse_matrix_for_each_triplet(rowIdx, transposeMatrixView, [&](IndexType row, IndexType col, ValueType value) {
            triplets.emplace_back(row, col, weightedValue(col, value) / rowSum);
        });
    }

    setMatrix(linalg::make_sparse_matrix_storage_from_triplets(transposeRows, transposeCols, triplets));
}


void Binning::print(std::ostream&) const {
    ATLAS_NOTIMPLEMENTED;
}

Binning::SparseMatrixStorage Binning::transposeAndHaloExchange(const SparseMatrixView& interpMatrix) const {
    // Transpose interpolation matrix. If transposed matrix contains rows which map into ghost elements,
    // move the rows to the PE which owns the element.

    struct Views {
        Views(const FunctionSpace& functionSpace):
            partition{array::make_view<int, 1>(functionSpace.partition())},
            remote_index{array::make_indexview<idx_t, 1>(functionSpace.remote_index())},
            ghost{array::make_view<int, 1>(functionSpace.ghost())} {}
        array::ArrayView<int, 1> partition;
        array::IndexView<idx_t, 1> remote_index;
        array::ArrayView<int, 1> ghost;
    };
    const Views sourceViews{source_};
    const Views targetViews{target_};

    // MPI send buffers.
    auto triplets   = std::vector<TripletType>{};
    auto sendBuffer = MpiBuffer{};

    linalg::sparse_matrix_for_each_triplet(interpMatrix, [&](IndexType row, IndexType col, ValueType weight) {
        // Swap row and column indices for transpose.
        std::swap(row, col);


        // Skip rows of interp matrix (i.e., cols of this transposed matrix) that are made redundant by halo exchanges.
        // These shouldn't exist, but some interpolation methods generate them anyway.
        if (sourceViews.ghost(col)) {
            return;
        }

        if (!targetViews.ghost(row)) {
            // Rows is owned by target function space.
            triplets.emplace_back(row, col, weight);
        }
        else {
            // Create MPI send data buffers for non-owned rows.
            const int remoteRank      = targetViews.partition(row);
            const IndexType remoteRow = static_cast<IndexType>(targetViews.remote_index(row));
            sendBuffer.pushBack(remoteRank, {remoteRow, col, weight});
        }
    });

    const auto recvBuffer = sendBuffer.allToAll();

    // Create a map from source partition and remote index to local index for ghost points only.
    auto sourceRemoteToLocalMap = std::map<std::pair<int, idx_t>, idx_t>{};
    for (idx_t localCol = 0; localCol < source_.size(); ++localCol) {
        if (sourceViews.ghost(localCol)) {
            const int remoteRank  = sourceViews.partition(localCol);
            const idx_t remoteCol = sourceViews.remote_index(localCol);
            sourceRemoteToLocalMap.insert({{remoteRank, remoteCol}, localCol});
        }
    }

    // Iterate over received data, creating owned rows for ghost columns
    for (std::size_t remoteRank = 0; remoteRank < mpi::comm().size(); ++remoteRank) {
        const auto subBuffer     = recvBuffer.get(remoteRank);
        const auto subBufferSize = recvBuffer.size(remoteRank);
        for (std::size_t i = 0; i < subBufferSize; ++i) {
            const auto triplet = subBuffer(i);

            const IndexType row       = triplet.row();
            const IndexType remoteCol = triplet.col();
            const ValueType weight    = triplet.value();

            ATLAS_ASSERT_MSG(!targetViews.ghost(row), "Row should be owned by the target functionspace");

            const auto col = [&] {
                try {
                    return sourceRemoteToLocalMap.at({remoteRank, remoteCol});
                }
                catch (const std::out_of_range&) {
                    ATLAS_THROW_EXCEPTION(
                        "Column local index not found. Try increasing source functionspace halo size.");
                }
            }();
            ATLAS_ASSERT_MSG(sourceViews.ghost(col), "Column should not be owned by the source functionspace.");
            triplets.emplace_back(row, col, weight);
        }
    }
    return linalg::make_sparse_matrix_storage_from_triplets(
        static_cast<IndexType>(target_.size()), static_cast<IndexType>(source_.size()), triplets);
}

std::vector<double> Binning::getAreaWeights() const {
    const auto csGrid    = CubedSphereGrid(source_.grid());
    auto ncFunctionSpace = functionspace::NodeColumns(source_);

    if (!(csGrid && ncFunctionSpace)) {
        // If not a cubed sphere grid and a node columns function space, return equal weights.
        return std::vector<double>(source_.size(), 1.);
    }

    auto mesh                  = ncFunctionSpace.mesh();
    const auto nodeWeights     = mesh::actions::GetCubedSphereNodalArea()(mesh);
    const auto nodeWeightsView = array::make_view<double, 1>(nodeWeights);

    auto weightVector = std::vector<double>(source_.size());
    for (idx_t i = 0; i < source_.size(); ++i) {
        weightVector[i] = nodeWeightsView(i);
    }
    return weightVector;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
