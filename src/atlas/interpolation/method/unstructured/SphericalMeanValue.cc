/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include <cmath>
#include <iomanip>
#include <limits>

#include "SphericalMeanValue.h"

#include "eckit/log/Plural.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/log/Seconds.h"

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"
#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/interpolation/element/Triag3D.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/Ray.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/Buffer.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"


namespace atlas {
namespace interpolation {
namespace method {

namespace {
MethodBuilder<SphericalMeanValue> __builder("spherical-mean-value");
}  // namespace

void SphericalMeanValue::do_setup(const Grid& source, const Grid& target, const Cache& cache) {
    allow_halo_exchange_ = false;
    //  no halo_exchange because we don't have any halo with delaunay or 3d structured meshgenerator

    if (interpolation::MatrixCache(cache)) {
        setMatrix(cache);
        ATLAS_ASSERT(matrix().rows() == target.size());
        ATLAS_ASSERT(matrix().cols() == source.size());
        return;
    }
    if (mpi::size() > 1) {
        ATLAS_NOTIMPLEMENTED;
    }
    auto make_nodecolumns = [](const Grid& grid) {
        Mesh mesh;
        if (StructuredGrid{grid}) {
            mesh = MeshGenerator("structured", util::Config("3d", true)).generate(grid);
        }
        else {
            mesh = MeshGenerator("delaunay").generate(grid);
        }
        return functionspace::NodeColumns(mesh);
    };

    do_setup(make_nodecolumns(source), functionspace::PointCloud{target});
}

void SphericalMeanValue::do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache& cache) {
    if (interpolation::MatrixCache(cache)) {
        setMatrix(cache);
        source_ = source;
        target_ = target;
        ATLAS_ASSERT(matrix().rows() == target.size());
        ATLAS_ASSERT(matrix().cols() == source.size());
        return;
    }
    do_setup(source, target);
}

void SphericalMeanValue::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    ATLAS_TRACE("atlas::interpolation::method::FiniteElement::do_setup()");

    source_ = source;
    target_ = target;

    ATLAS_TRACE_SCOPE("Setup target") {
        auto create_xyz = [](Field lonlat_field) {
            auto xyz_field = Field("xyz", array::make_datatype<double>(), array::make_shape(lonlat_field.shape(0), 3));
            auto lonlat    = array::make_view<double, 2>(lonlat_field);
            auto xyz       = array::make_view<double, 2>(xyz_field);
            PointXYZ p2;
            for (idx_t n = 0; n < lonlat.shape(0); ++n) {
                const PointLonLat p1(lonlat(n, 0), lonlat(n, 1));
                util::Earth::convertSphericalToCartesian(p1, p2);
                xyz(n, 0) = p2.x();
                xyz(n, 1) = p2.y();
                xyz(n, 2) = p2.z();
            }
            return xyz_field;
        };

        target_ghost_  = target.ghost();
        target_lonlat_ = target.lonlat();
        if (functionspace::NodeColumns tgt = target) {
            auto meshTarget = tgt.mesh();
            target_xyz_     = mesh::actions::BuildXYZField("xyz")(meshTarget);
        }
        else {
            target_xyz_ = create_xyz(target_lonlat_);
        }
    }

    setup(source);
}

struct Stencil {
    enum
    {
        max_stencil_size = 4
    };
};

void SphericalMeanValue::print(std::ostream& out) const {
    functionspace::NodeColumns src(source_);
    functionspace::NodeColumns tgt(target_);
    out << "atlas::interpolation::method::FiniteElement{" << std::endl;
    out << "max_fraction_elems_to_try: " << max_fraction_elems_to_try_;
    out << ", treat_failure_as_missing_value: " << treat_failure_as_missing_value_;
    if (not tgt) {
        out << "}" << std::endl;
        return;
    }
    out << ", NodeColumns to NodeColumns stencil weights: " << std::endl;
    auto gidx_src = array::make_view<gidx_t, 1>(src.nodes().global_index());

    ATLAS_ASSERT(tgt.nodes().size() == idx_t(matrix().rows()));


    auto field_stencil_points_loc  = tgt.createField<gidx_t>(option::variables(Stencil::max_stencil_size));
    auto field_stencil_weights_loc = tgt.createField<double>(option::variables(Stencil::max_stencil_size));
    auto field_stencil_size_loc    = tgt.createField<int>();

    auto stencil_points_loc  = array::make_view<gidx_t, 2>(field_stencil_points_loc);
    auto stencil_weights_loc = array::make_view<double, 2>(field_stencil_weights_loc);
    auto stencil_size_loc    = array::make_view<idx_t, 1>(field_stencil_size_loc);
    stencil_size_loc.assign(0);

    const auto m = atlas::linalg::make_non_owning_eckit_sparse_matrix(matrix());
    for (auto it = m.begin(); it != m.end(); ++it) {
        idx_t p                   = idx_t(it.row());
        idx_t& i                  = stencil_size_loc(p);
        stencil_points_loc(p, i)  = gidx_src(it.col());
        stencil_weights_loc(p, i) = *it;
        ++i;
    }

    gidx_t global_size = tgt.gather().glb_dof();

    auto field_stencil_points_glb =
        tgt.createField<gidx_t>(option::variables(Stencil::max_stencil_size) | option::global(0));
    auto field_stencil_weights_glb =
        tgt.createField<double>(option::variables(Stencil::max_stencil_size) | option::global(0));
    auto field_stencil_size_glb = tgt.createField<idx_t>(option::global(0));


    auto stencil_points_glb  = array::make_view<gidx_t, 2>(field_stencil_points_glb);
    auto stencil_weights_glb = array::make_view<double, 2>(field_stencil_weights_glb);
    auto stencil_size_glb    = array::make_view<idx_t, 1>(field_stencil_size_glb);

    tgt.gather().gather(stencil_size_loc, stencil_size_glb);
    tgt.gather().gather(stencil_points_loc, stencil_points_glb);
    tgt.gather().gather(stencil_weights_loc, stencil_weights_glb);

    int precision = std::numeric_limits<double>::max_digits10;
    for (idx_t i = 0; i < global_size; ++i) {
        out << std::setw(10) << i + 1 << " : ";
        for (idx_t j = 0; j < stencil_size_glb(i); ++j) {
            out << std::setw(10) << stencil_points_glb(i, j);
        }
        for (idx_t j = stencil_size_glb(i); j < Stencil::max_stencil_size; ++j) {
            out << "          ";
        }
        for (idx_t j = 0; j < stencil_size_glb(i); ++j) {
            out << std::setw(precision + 5) << std::left << std::setprecision(precision) << stencil_weights_glb(i, j);
        }
        out << std::endl;
    }
    out << "}" << std::endl;
}

void SphericalMeanValue::setup(const FunctionSpace& source) {
    const functionspace::NodeColumns src = source;
    ATLAS_ASSERT(src);

    ocoords_.reset(new array::ArrayView<double, 2>(array::make_view<double, 2>(target_xyz_)));
    idx_t out_npts = ocoords_->shape(0);
    // return early if no output points on this partition reserve is called on
    // the triplets but also during the sparseMatrix constructor. This won't
    // work for empty matrices
    if (out_npts == 0) {
        return;
    }

    Mesh meshSource = src.mesh();


    auto trace_setup_source = atlas::Trace{Here(), "Setup source"};

    // generate 3D point coordinates
    Field source_xyz = mesh::actions::BuildXYZField("xyz")(meshSource);

    // generate barycenters of each triangle & insert them on a kd-tree
    util::Config config;
    config.set("name", "centre ");
    config.set("flatten_virtual_elements", false);
    Field cell_centres = mesh::actions::BuildCellCentres(config)(meshSource);

    std::unique_ptr<ElemIndex3> eTree(create_element_kdtree(meshSource, cell_centres));

    trace_setup_source.stop();


    icoords_.reset(new array::ArrayView<double, 2>(array::make_view<double, 2>(source_xyz)));
    igidx_.reset(new array::ArrayView<gidx_t, 1>(array::make_view<gidx_t, 1>(src.nodes().global_index())));
    connectivity_              = &meshSource.cells().node_connectivity();
    const mesh::Nodes& i_nodes = meshSource.nodes();

    idx_t inp_npts = i_nodes.size();

    auto target_point = [this](idx_t ip) {
        return PointXYZ{(*ocoords_)(ip, 0), (*ocoords_)(ip, 1), (*ocoords_)(ip, 2)};
    };


    array::ArrayView<int, 1> out_ghosts = array::make_view<int, 1>(target_ghost_);

    array::ArrayView<double, 2> out_lonlat = array::make_view<double, 2>(target_lonlat_);

    idx_t Nelements = meshSource.cells().size();

    // weights -- one per vertex of element, triangles (3) or quads (4)

    Triplets weights_triplets;              // structure to fill-in sparse matrix
    weights_triplets.resize(out_npts * 4);  // preallocate space as if all elements where quads
    auto insert_triplets = [&weights_triplets](idx_t n, const Triplets& triplets) -> bool {
        if (triplets.size()) {
            std::copy(triplets.begin(), triplets.end(), weights_triplets.begin() + 4 * n);
            return true;
        }
        return false;
    };

    double search_radius = 0.;
    if (meshSource.metadata().has("cell_maximum_diagonal_on_unit_sphere")) {
        search_radius =
            Geometry("Earth").radius() * meshSource.metadata().getDouble("cell_maximum_diagonal_on_unit_sphere");
        ASSERT(search_radius > 0.);
        Log::debug() << "k-d tree: search radius = " << search_radius / 1000. << " km" << std::endl;
    }
    auto find_element_candidates_in_search_radius = [&eTree, &search_radius](const PointXYZ& p) {
        return eTree->findInSphere(p, search_radius);
    };
    auto find_k_nearest_element_candidates = [&eTree](const PointXYZ& p, size_t k) {
        return eTree->kNearestNeighbours(p, k);
    };
    auto try_interpolate_with_element_candidates =
        [&insert_triplets, this](idx_t n, const ElemIndex3::NodeList& element_candidates) -> bool {
        if (element_candidates.empty()) {
            return false;
        }
        return insert_triplets(n, projectPointToElements(n, element_candidates));
    };


    // search nearest k cell centres
    const idx_t maxNbElemsToTry             = std::max<idx_t>(8, idx_t(Nelements * max_fraction_elems_to_try_));
    size_t diagnosed_max_neighbours         = 0;
    bool allowed_to_diagnose_max_neighbours = Log::debug() && atlas_omp_get_max_threads() > 1;

    std::vector<size_t> failures;

    ATLAS_TRACE_SCOPE("Computing interpolation matrix") {
        std::unique_ptr<eckit::ProgressTimer> progress;
        if (atlas_omp_get_max_threads() == 1) {
            progress.reset(new eckit::ProgressTimer{"Computing interpolation weights", static_cast<size_t>(out_npts),
                                                    "point", double(5), Log::debug()});
        }
        atlas_omp_parallel_for(idx_t ip = 0; ip < out_npts; ++ip) {
            if (out_ghosts(ip)) {
                continue;
            }

            bool success = false;
            if (search_radius != 0.) {
                auto p  = target_point(ip);
                success = try_interpolate_with_element_candidates(ip, find_element_candidates_in_search_radius(p));
            }
            else {
                size_t k = 1;
                auto p   = target_point(ip);
                while (!success && k <= maxNbElemsToTry) {
                    if (allowed_to_diagnose_max_neighbours) {  // avoid race condition
                        diagnosed_max_neighbours = std::max(k, diagnosed_max_neighbours);
                    }
                    success = try_interpolate_with_element_candidates(ip, find_k_nearest_element_candidates(p, k));
                    k *= 2;
                }
            }

            if (!success) {
                atlas_omp_critical {
                    failures.emplace_back(ip);
                }
            }
            if (progress) {
                ++(*progress);
            }
        }
    }
    if (diagnosed_max_neighbours) {
        Log::debug() << "Maximum neighbours searched was " << eckit::Plural(diagnosed_max_neighbours, "element")
                     << std::endl;
    }

    if (failures.size()) {
        if (treat_failure_as_missing_value_) {
            missing_.resize(failures.size());
            std::copy(std::begin(failures), std::end(failures), missing_.begin());
        }
        else {
            // If this fails, consider lowering atlas::grid::parametricEpsilon
            std::ostringstream msg;
            msg << "Rank " << mpi::rank() << " failed to project points:\n";
            for (std::vector<size_t>::const_iterator i = failures.begin(); i != failures.end(); ++i) {
                const PointLonLat pll{out_lonlat(*i, (size_t)0), out_lonlat(*i, (size_t)1)};  // lookup point
                msg << "\t(lon,lat) = " << pll << "\n";
            }

            Log::error() << msg.str() << std::endl;
            throw_Exception(msg.str());
        }
    }

    // fill sparse matrix and return, this cannot be multithreaded!
    setMatrix(out_npts, inp_npts, weights_triplets);
}

Method::Triplets SphericalMeanValue::projectPointToElements(size_t ip, const ElemIndex3::NodeList& elems) const {
    ATLAS_ASSERT(elems.begin() != elems.end());

    const size_t inp_points = icoords_->shape(0);

    Triplets triplets;
    triplets.reserve(4);

    const PointXYZ candidatePoint = PointXYZ::normalize(
        PointXYZ{(*ocoords_)(ip, size_t(0)), (*ocoords_)(ip, size_t(1)), (*ocoords_)(ip, size_t(2))});

    for (ElemIndex3::NodeList::const_iterator itc = elems.begin(); itc != elems.end(); ++itc) {
        const idx_t elem_id = idx_t((*itc).value().payload());
        ATLAS_ASSERT(elem_id < connectivity_->rows());

        size_t nb_cols = connectivity_->cols(elem_id);

        std::vector<size_t> idx;
        idx.reserve(nb_cols);

        for (idx_t i = 0; i < nb_cols; ++i) {
            idx.emplace_back((*connectivity_)(elem_id, i));
            ATLAS_ASSERT(idx[i] < inp_points);
        }

        std::vector<PointXYZ> listVertices;
        listVertices.reserve(util::ConvexSphericalPolygon::MAX_SIZE);

        for (idx_t i = 0; i < nb_cols; ++i) {
            listVertices.emplace_back(PointXYZ::normalize(PointXYZ{
                (*icoords_)(idx[i], size_t(0)), (*icoords_)(idx[i], size_t(1)), (*icoords_)(idx[i], size_t(2))}));
        }

        util::ConvexSphericalPolygon currentPolygon(listVertices.data(), listVertices.size());

        nb_cols = currentPolygon.size();
        ATLAS_ASSERT(nb_cols == 3 || nb_cols == 4);

        std::vector<double> polygonWeights(nb_cols);

        if (currentPolygon.compute_vertex_weights(candidatePoint, polygonWeights.data(), polygonWeights.size()) == 1) {
            if (normalisation_) {
                double weightsSum = 0;
                for (size_t i = 0; i < nb_cols; ++i) {
                    weightsSum += (polygonWeights)[i];
                }
                for (size_t i = 0; i < nb_cols; ++i) {
                    (polygonWeights)[i] /= weightsSum;
                }
            }
            for (size_t i = 0; i < nb_cols; ++i) {
                triplets.emplace_back(ip, idx[i], (polygonWeights)[i]);
            }
            break;
        }
    }

    return triplets;
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
