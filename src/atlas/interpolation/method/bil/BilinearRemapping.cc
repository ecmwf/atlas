/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iomanip>
#include <limits>

#include "BilinearRemapping.h"

#include "eckit/log/Plural.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/log/Seconds.h"

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"
#include "atlas/interpolation/element/Quad2D.h"
#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/Ray.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/Build2DCellCentres.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/Buffer.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"


namespace atlas {
namespace interpolation {
namespace method {

namespace {

MethodBuilder<BilinearRemapping> __builder("bilinear-remapping");

// epsilon used to scale edge tolerance when projecting ray to intesect element
static const double parametricEpsilon = 1e-15;

}  // namespace


void BilinearRemapping::do_setup(const Grid& source, const Grid& target, const Cache& cache) {
    allow_halo_exchange_ = false;
    //  no halo_exchange because we don't have any halo with delaunay or 3d structured meshgenerator

    if (interpolation::MatrixCache(cache)) {
        matrix_cache_ = cache;
        matrix_       = &matrix_cache_.matrix();
        ATLAS_ASSERT(matrix_cache_.matrix().rows() == target.size());
        ATLAS_ASSERT(matrix_cache_.matrix().cols() == source.size());
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

void BilinearRemapping::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    ATLAS_TRACE("atlas::interpolation::method::BilinearRemapping::do_setup()");

    source_ = source;
    target_ = target;

    ATLAS_TRACE_SCOPE("Setup target") {
        if (functionspace::NodeColumns tgt = target) {
            Mesh meshTarget = tgt.mesh();

            target_ghost_  = meshTarget.nodes().ghost();
            target_lonlat_ = meshTarget.nodes().lonlat();
        }
        else if (functionspace::PointCloud tgt = target) {
            const idx_t N  = tgt.size();
            target_ghost_  = tgt.ghost();
            target_lonlat_ = tgt.lonlat();
        }
        else {
            ATLAS_NOTIMPLEMENTED;
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

void BilinearRemapping::print(std::ostream& out) const {
    functionspace::NodeColumns src(source_);
    functionspace::NodeColumns tgt(target_);
    out << "atlas::interpolation::method::BilinearRemapping{" << std::endl;
    out << "max_fraction_elems_to_try: " << max_fraction_elems_to_try_;
    out << ", treat_failure_as_missing_value: " << treat_failure_as_missing_value_;
    if (not tgt) {
        out << "}" << std::endl;
        return;
    }
    out << ", NodeColumns to NodeColumns stencil weights: " << std::endl;
    auto gidx_src = array::make_view<gidx_t, 1>(src.nodes().global_index());

    ATLAS_ASSERT(tgt.nodes().size() == idx_t(matrix_->rows()));


    auto field_stencil_points_loc  = tgt.createField<gidx_t>(option::variables(Stencil::max_stencil_size));
    auto field_stencil_weights_loc = tgt.createField<double>(option::variables(Stencil::max_stencil_size));
    auto field_stencil_size_loc    = tgt.createField<int>();

    auto stencil_points_loc  = array::make_view<gidx_t, 2>(field_stencil_points_loc);
    auto stencil_weights_loc = array::make_view<double, 2>(field_stencil_weights_loc);
    auto stencil_size_loc    = array::make_view<idx_t, 1>(field_stencil_size_loc);
    stencil_size_loc.assign(0);

    for (Matrix::const_iterator it = matrix_->begin(); it != matrix_->end(); ++it) {
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

void BilinearRemapping::setup(const FunctionSpace& source) {
    const functionspace::NodeColumns src = source;
    ATLAS_ASSERT(src);

    Mesh meshSource = src.mesh();


    auto trace_setup_source = atlas::Trace{Here(), "Setup source"};

    // 2D point coordinates
    auto source_lonlat = array::make_view<double, 2>(src.lonlat());

    // generate barycenters of each triangle & insert them on a kd-tree
    util::Config config;
    config.set("name", "centre ");
    config.set("flatten_virtual_elements", false);
    Field cell_centres = mesh::actions::Build2DCellCentres(config)(meshSource);

    std::unique_ptr<ElemIndex2> eTree(create_element2D_kdtree(meshSource, cell_centres));

    trace_setup_source.stop();

    ilonlat_.reset(new array::ArrayView<double, 2>(array::make_view<double, 2>(meshSource.nodes().lonlat())));
    olonlat_.reset(new array::ArrayView<double, 2>(array::make_view<double, 2>(target_lonlat_)));
    igidx_.reset(new array::ArrayView<gidx_t, 1>(array::make_view<gidx_t, 1>(src.nodes().global_index())));
    connectivity_              = &meshSource.cells().node_connectivity();
    const mesh::Nodes& i_nodes = meshSource.nodes();


    idx_t inp_npts = i_nodes.size();
    idx_t out_npts = olonlat_->shape(0);

    array::ArrayView<int, 1> out_ghosts = array::make_view<int, 1>(target_ghost_);

    idx_t Nelements = meshSource.cells().size();

    // weights -- one per vertex of element, triangles (3) or quads (4)

    Triplets weights_triplets;               // structure to fill-in sparse matrix
    weights_triplets.reserve(out_npts * 4);  // preallocate space as if all elements where quads

    // search nearest k cell centres

    const idx_t maxNbElemsToTry = std::max<idx_t>(8, idx_t(Nelements * max_fraction_elems_to_try_));
    idx_t max_neighbours        = 0;

    std::vector<size_t> failures;

    ATLAS_TRACE_SCOPE("Computing interpolation matrix") {
        eckit::ProgressTimer progress("Computing interpolation weights", out_npts, "point", double(5), Log::debug());
        for (idx_t ip = 0; ip < out_npts; ++ip, ++progress) {
            if (out_ghosts(ip)) {
                continue;
            }

            double lon{(*olonlat_)(ip, 0)};
            while (lon > 360.0) {
                lon -= 360.0;
            }
            while (lon < 0.0) {
                lon += 360.0;
            }
            PointXY p{lon, (*olonlat_)(ip, 1)};  // lookup point

            idx_t kpts   = 1;
            bool success = false;
            std::ostringstream failures_log;

            while (!success && kpts <= maxNbElemsToTry) {
                max_neighbours = std::max(kpts, max_neighbours);

                ElemIndex2::NodeList cs = eTree->kNearestNeighbours(p, kpts);
                Triplets triplets       = projectPointToElements(ip, cs, failures_log);

                if (triplets.size()) {
                    std::copy(triplets.begin(), triplets.end(), std::back_inserter(weights_triplets));
                    success = true;
                }
                kpts *= 2;
            }

            if (!success) {
                failures.push_back(ip);
                if (not treat_failure_as_missing_value_) {
                    Log::debug() << "------------------------------------------------------"
                                    "---------------------\n";
                    Log::debug() << "Failed to project point (lon,lat)=" << p << '\n';
                    Log::debug() << failures_log.str();
                }
            }
        }
    }
    Log::debug() << "Maximum neighbours searched was " << eckit::Plural(max_neighbours, "element") << std::endl;

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
                const PointLonLat pll{(*olonlat_)(*i, (size_t)0), (*olonlat_)(*i, (size_t)1)};  // lookup point
                msg << "\t(lon,lat) = " << pll << "\n";
            }

            Log::error() << msg.str() << std::endl;
            throw_Exception(msg.str());
        }
    }

    // fill sparse matrix and return
    Matrix A(out_npts, inp_npts, weights_triplets);
    matrix_shared_->swap(A);
}

Method::Triplets BilinearRemapping::projectPointToElements(size_t ip, const ElemIndex2::NodeList& elems,
                                                           std::ostream& /* failures_log */) const {
    ATLAS_ASSERT(elems.begin() != elems.end());

    const size_t inp_points = ilonlat_->shape(0);
    std::array<size_t, 4> idx;
    std::array<double, 4> w;

    Triplets triplets;
    triplets.reserve(4);

    double lon{(*olonlat_)(ip, 0)};
    while (lon > 360.0) {
        lon -= 360.0;
    }
    while (lon < 0.0) {
        lon += 360.0;
    }
    PointXY ob_loc{lon, (*olonlat_)(ip, 1)};  // lookup point

    for (ElemIndex2::NodeList::const_iterator itc = elems.begin(); itc != elems.end(); ++itc) {
        const idx_t elem_id = idx_t((*itc).value().payload());
        ATLAS_ASSERT(elem_id < connectivity_->rows());

        const idx_t nb_cols = connectivity_->cols(elem_id);
        ATLAS_ASSERT(nb_cols == 3 || nb_cols == 4);

        for (idx_t i = 0; i < nb_cols; ++i) {
            idx[i] = (*connectivity_)(elem_id, i);
            ATLAS_ASSERT(idx[i] < inp_points);
        }

        if (nb_cols == 3) {
            /* triangle */
            element::Triag2D triag(PointXY{(*ilonlat_)(idx[0], 0), (*ilonlat_)(idx[0], 1)},
                                   PointXY{(*ilonlat_)(idx[1], 0), (*ilonlat_)(idx[1], 1)},
                                   PointXY{(*ilonlat_)(idx[2], 0), (*ilonlat_)(idx[2], 1)});

            // pick an epsilon based on a characteristic length (sqrt(area))
            // (this scales linearly so it better compares with linear weights u,v,w)
            const double edgeEpsilon = parametricEpsilon * std::sqrt(triag.area());
            ATLAS_ASSERT(edgeEpsilon >= 0);

            Intersect is = triag.intersects(ob_loc, edgeEpsilon);

            if (is) {
                // weights are the linear Lagrange function evaluated at u,v (aka
                // barycentric coordinates)
                w[0] = 1. - is.u - is.v;
                w[1] = is.u;
                w[2] = is.v;

                for (size_t i = 0; i < 3; ++i) {
                    triplets.emplace_back(ip, idx[i], w[i]);
                }

                break;  // stop looking for elements
            }
        }
        else {
            /* quadrilateral */
            element::Quad2D quad(PointXY{(*ilonlat_)(idx[0], 0), (*ilonlat_)(idx[0], 1)},
                                 PointXY{(*ilonlat_)(idx[1], 0), (*ilonlat_)(idx[1], 1)},
                                 PointXY{(*ilonlat_)(idx[2], 0), (*ilonlat_)(idx[2], 1)},
                                 PointXY{(*ilonlat_)(idx[3], 0), (*ilonlat_)(idx[3], 1)});

            // pick an epsilon based on a characteristic length (sqrt(area))
            // (this scales linearly so it better compares with linear weights u,v,w)
            const double edgeEpsilon = parametricEpsilon * std::sqrt(quad.area());
            ATLAS_ASSERT(edgeEpsilon >= 0);

            Intersect is = quad.localRemap(ob_loc, edgeEpsilon);

            if (is) {
                //std::cout << "intersection found for p: ";
                //std::cout << ob_loc << " and " << quad << std::endl;
                // weights are the bilinear Lagrange function evaluated at u,v
                w[0] = (1. - is.u) * (1. - is.v);
                w[1] = is.u * (1. - is.v);
                w[2] = is.u * is.v;
                w[3] = (1. - is.u) * is.v;

                for (size_t i = 0; i < 4; ++i) {
                    triplets.emplace_back(ip, idx[i], w[i]);
                }
                break;  // stop looking for elements
            }
        }

    }  // loop over nearest elements

    if (!triplets.empty()) {
        normalise(triplets);
    }
    return triplets;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
