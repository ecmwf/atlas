/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iomanip>

#include "eckit/log/Bytes.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/detail/MeshImpl.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Unique.h"

namespace atlas {
namespace mesh {
namespace detail {

//----------------------------------------------------------------------------------------------------------------------

PartitionGraph* build_partition_graph(const MeshImpl& mesh) {
    ATLAS_TRACE("build_partition_graph...");
    const eckit::mpi::Comm& comm = mpi::comm();
    const int mpi_size           = int(comm.size());

    const util::Polygon& poly = mesh.polygon();

    std::vector<double> polygon;
    polygon.reserve(poly.size() * 2);

    auto xy = array::make_view<double, 2>(mesh.nodes().xy());

    for (idx_t node : poly) {
        polygon.push_back(xy(node, XX));
        polygon.push_back(xy(node, YY));
    }

    eckit::mpi::Buffer<double> recv_polygons(mpi_size);
    comm.allGatherv(polygon.begin(), polygon.end(), recv_polygons);

    using PolygonXY = std::vector<PointXY>;
    std::vector<PolygonXY> polygons(mpi_size);
    for (idx_t p = 0; p < mpi_size; ++p) {
        for (idx_t j = 0; j < recv_polygons.counts[p] / 2; ++j) {
            PointXY pxy(*(recv_polygons.begin() + recv_polygons.displs[p] + 2 * j + XX),
                        *(recv_polygons.begin() + recv_polygons.displs[p] + 2 * j + YY));
            polygons[p].push_back(pxy);
        }
    }

    std::map<uidx_t, std::set<idx_t>> uid_2_parts;
    idx_t jpart = 0;
    for (const PolygonXY& _polygon : polygons) {
        for (const PointXY& pxy : _polygon) {
            PointLonLat pll = pxy;
            if (eckit::types::is_strictly_greater(0., pll.lon())) {
                pll.lon() += 360.;
            }
            if (eckit::types::is_approximately_greater_or_equal(pll.lon(), 360.)) {
                pll.lon() -= 360.;
            }
            uidx_t uid = util::unique_lonlat(pll.data());
            uid_2_parts[uid].insert(jpart);
        }
        ++jpart;
    }
    std::vector<std::set<idx_t>> graph(mpi_size);
    for (const auto& u2p : uid_2_parts) {
        const std::set<idx_t>& parts = u2p.second;
        for (idx_t jpart : parts) {
            for (idx_t ipart : parts) {
                if (jpart != ipart) {
                    graph[jpart].insert(ipart);
                }
            }
        }
    }

    std::vector<idx_t> counts(mpi_size);
    std::vector<idx_t> displs(mpi_size);
    idx_t values_size = 0;
    for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
        counts[jpart] = graph[jpart].size();
        displs[jpart] = values_size;
        values_size += counts[jpart];
    }
    std::vector<idx_t> values;
    values.reserve(values_size);
    for (const std::set<idx_t>& graph_node : graph) {
        for (idx_t v : graph_node) {
            values.push_back(v);
        }
    }

    return new PartitionGraph(values.data(), mpi_size, displs.data(), counts.data());
}

size_t PartitionGraph::footprint() const {
    size_t size = sizeof(*this);
    size += sizeof(idx_t) * displs_.capacity();
    size += sizeof(idx_t) * counts_.capacity();
    size += sizeof(idx_t) * values_.capacity();
    return size;
}

idx_t PartitionGraph::size() const {
    return displs_.size();
}

PartitionGraph::Neighbours PartitionGraph::nearestNeighbours(const idx_t partition) const {
    return Neighbours(values_.data() + displs_[partition], values_.data() + displs_[partition] + counts_[partition]);
}

PartitionGraph::PartitionGraph() = default;

PartitionGraph::PartitionGraph(idx_t values[], idx_t rows, idx_t displs[], idx_t counts[]) {
    displs_.assign(displs, displs + rows);
    counts_.assign(counts, counts + rows);
    values_.assign(values, values + displs[rows - 1] + counts[rows - 1]);

    for (idx_t jpart = 0; jpart < rows; ++jpart) {
        for (idx_t neighbour : nearestNeighbours(jpart)) {
            bool found(false);
            for (idx_t nextneighbour : nearestNeighbours(neighbour)) {
                if (nextneighbour == jpart) {
                    found = true;
                }
            }
            if (not found) {
                values_.insert(values_.begin() + displs_[neighbour] + counts_[neighbour], jpart);
                counts_[neighbour]++;
                for (idx_t j = neighbour + 1; j < rows; ++j) {
                    displs_[j]++;
                }
            }
        }
    }

    maximum_nearest_neighbours_ = 0;
    for (idx_t n : counts_) {
        maximum_nearest_neighbours_ = std::max(n, maximum_nearest_neighbours_);
    }
}

idx_t PartitionGraph::maximumNearestNeighbours() const {
    return maximum_nearest_neighbours_;
}

void PartitionGraph::print(std::ostream& os) const {
    for (idx_t jpart = 0; jpart < size(); ++jpart) {
        os << std::setw(3) << jpart << " : ";
        for (idx_t v : nearestNeighbours(jpart)) {
            os << std::setw(3) << v << " ";
        }
        os << '\n';
    }
    os << "partition graph maximum neighbours = " << maximumNearestNeighbours() << '\n';
    os << "partition graph footprint = " << eckit::Bytes(footprint());
}

PartitionGraph::operator bool() const {
    return size();
}

std::ostream& operator<<(std::ostream& s, const PartitionGraph& p) {
    p.print(s);
    return s;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace mesh
}  // namespace atlas
