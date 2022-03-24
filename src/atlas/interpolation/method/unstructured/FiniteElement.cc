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

#include "FiniteElement.h"

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
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildXYZField.h"
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

MethodBuilder<FiniteElement> __builder("finite-element");

// epsilon used to scale edge tolerance when projecting ray to intesect element
static const double parametricEpsilon = 1e-15;

}  // namespace


void FiniteElement::do_setup(const Grid& source, const Grid& target, const Cache& cache) {
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

void FiniteElement::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    ATLAS_TRACE("atlas::interpolation::method::FiniteElement::do_setup()");

    source_ = source;
    target_ = target;

    ATLAS_TRACE_SCOPE("Setup target") {
        if (functionspace::NodeColumns tgt = target) {
            Mesh meshTarget = tgt.mesh();

            // generate 3D point coordinates
            target_xyz_    = mesh::actions::BuildXYZField("xyz")(meshTarget);
            target_ghost_  = meshTarget.nodes().ghost();
            target_lonlat_ = meshTarget.nodes().lonlat();
        }
        else if (functionspace::PointCloud tgt = target) {
            const idx_t N  = tgt.size();
            target_xyz_    = Field("xyz", array::make_datatype<double>(), array::make_shape(N, 3));
            target_ghost_  = tgt.ghost();
            target_lonlat_ = tgt.lonlat();
            auto lonlat    = array::make_view<double, 2>(tgt.lonlat());
            auto xyz       = array::make_view<double, 2>(target_xyz_);
            PointXYZ p2;
            for (idx_t n = 0; n < N; ++n) {
                const PointLonLat p1(lonlat(n, 0), lonlat(n, 1));
                util::Earth::convertSphericalToCartesian(p1, p2);
                xyz(n, 0) = p2.x();
                xyz(n, 1) = p2.y();
                xyz(n, 2) = p2.z();
            }
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

void FiniteElement::print(std::ostream& out) const {
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

    for (auto it = matrix().begin(); it != matrix().end(); ++it) {
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

void FiniteElement::setup(const FunctionSpace& source) {
    const functionspace::NodeColumns src = source;
    ATLAS_ASSERT(src);

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
    ocoords_.reset(new array::ArrayView<double, 2>(array::make_view<double, 2>(target_xyz_)));
    igidx_.reset(new array::ArrayView<gidx_t, 1>(array::make_view<gidx_t, 1>(src.nodes().global_index())));
    connectivity_              = &meshSource.cells().node_connectivity();
    const mesh::Nodes& i_nodes = meshSource.nodes();


    idx_t inp_npts = i_nodes.size();
    idx_t out_npts = ocoords_->shape(0);

    array::ArrayView<int, 1> out_ghosts = array::make_view<int, 1>(target_ghost_);

    array::ArrayView<double, 2> out_lonlat = array::make_view<double, 2>(target_lonlat_);

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

            PointXYZ p{(*ocoords_)(ip, 0), (*ocoords_)(ip, 1), (*ocoords_)(ip, 2)};  // lookup point

            idx_t kpts   = 1;
            bool success = false;
            std::ostringstream failures_log;

            while (!success && kpts <= maxNbElemsToTry) {
                max_neighbours = std::max(kpts, max_neighbours);

                ElemIndex3::NodeList cs = eTree->kNearestNeighbours(p, kpts);
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
                    const PointLonLat pll{out_lonlat(ip, 0), out_lonlat(ip, 1)};
                    Log::debug() << "Failed to project point (lon,lat)=" << pll << '\n';
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
                const PointLonLat pll{out_lonlat(*i, (size_t)0), out_lonlat(*i, (size_t)1)};  // lookup point
                msg << "\t(lon,lat) = " << pll << "\n";
            }

            Log::error() << msg.str() << std::endl;
            throw_Exception(msg.str());
        }
    }

    // fill sparse matrix and return
    Matrix A(out_npts, inp_npts, weights_triplets);
    setMatrix(A);
}

struct ElementEdge {
    std::array<idx_t, 2> idx;
    void swap() {
        idx_t tmp = idx[0];
        idx[0]    = idx[1];
        idx[1]    = tmp;
    }
};

Method::Triplets FiniteElement::projectPointToElements(size_t ip, const ElemIndex3::NodeList& elems,
                                                       std::ostream& /* failures_log */) const {
    ATLAS_ASSERT(elems.begin() != elems.end());

    const size_t inp_points = icoords_->shape(0);
    std::array<size_t, 4> idx;
    std::array<double, 4> w;

    Triplets triplets;
    triplets.reserve(4);
    Ray ray(PointXYZ{(*ocoords_)(ip, size_t(0)), (*ocoords_)(ip, size_t(1)), (*ocoords_)(ip, size_t(2))});
    const Vector3D p{(*ocoords_)(ip, size_t(0)), (*ocoords_)(ip, size_t(1)), (*ocoords_)(ip, size_t(2))};
    ElementEdge edge;
    idx_t single_point;
    for (ElemIndex3::NodeList::const_iterator itc = elems.begin(); itc != elems.end(); ++itc) {
        const idx_t elem_id = idx_t((*itc).value().payload());
        ATLAS_ASSERT(elem_id < connectivity_->rows());

        const idx_t nb_cols = connectivity_->cols(elem_id);
        ATLAS_ASSERT(nb_cols == 3 || nb_cols == 4);

        for (idx_t i = 0; i < nb_cols; ++i) {
            idx[i] = (*connectivity_)(elem_id, i);
            ATLAS_ASSERT(idx[i] < inp_points);
        }

        constexpr double tolerance = 1.e-12;

        auto on_triag_edge = [&]() {
            if (w[0] < tolerance) {
                edge.idx[0] = 1;
                edge.idx[1] = 2;
                w[0]        = 0.;
                return true;
            }
            if (w[1] < tolerance) {
                edge.idx[0] = 0;
                edge.idx[1] = 2;
                w[1]        = 0.;
                return true;
            }
            if (w[2] < tolerance) {
                edge.idx[0] = 0;
                edge.idx[1] = 1;
                w[2]        = 0.;
                return true;
            }
            return false;
        };

        auto on_quad_edge = [&]() {
            if (w[0] < tolerance && w[1] < tolerance) {
                edge.idx[0] = 2;
                edge.idx[1] = 3;
                w[0]        = 0.;
                w[1]        = 0.;
                return true;
            }
            if (w[1] < tolerance && w[2] < tolerance) {
                edge.idx[0] = 0;
                edge.idx[1] = 3;
                w[1]        = 0.;
                w[2]        = 0.;
                return true;
            }
            if (w[2] < tolerance && w[3] < tolerance) {
                edge.idx[0] = 0;
                edge.idx[1] = 1;
                w[2]        = 0.;
                w[3]        = 0.;
                return true;
            }
            if (w[3] < tolerance && w[0] < tolerance) {
                edge.idx[0] = 1;
                edge.idx[1] = 2;
                w[3]        = 0.;
                w[0]        = 0.;
                return true;
            }
            return false;
        };

        auto on_single_point = [&]() {
            if (w[edge.idx[0]] < tolerance) {
                single_point   = edge.idx[1];
                w[edge.idx[0]] = 0.;
                return true;
            }
            if (w[edge.idx[1]] < tolerance) {
                single_point   = edge.idx[0];
                w[edge.idx[1]] = 0.;
                return true;
            }
            return false;
        };

        auto interpolate_edge = [&](const Vector3D& p0, const Vector3D& p1) {
            /*
             * Given points p0,p1 defining the edge, and point p, find projected point pt
             * on edge to compute interpolation weights.
             *                  p
             *                  |`.
             *                  |  `.v
             *                  |    `.
             *  p1--------------pt-----p0
             *                  <--d----
             */
            Vector3D d     = (p1 - p0) / (p1 - p0).norm();
            Vector3D v     = p - p0;
            double t       = v.dot(d);
            Vector3D pt    = p0 + d * t;
            t              = (pt - p0).norm() / (p1 - p0).norm();
            w[edge.idx[0]] = 1. - t;
            w[edge.idx[1]] = t;
        };

        if (nb_cols == 3) {
            /* triangle */
            element::Triag3D triag(PointXYZ{(*icoords_)(idx[0], size_t(0)), (*icoords_)(idx[0], size_t(1)),
                                            (*icoords_)(idx[0], size_t(2))},
                                   PointXYZ{(*icoords_)(idx[1], size_t(0)), (*icoords_)(idx[1], size_t(1)),
                                            (*icoords_)(idx[1], size_t(2))},
                                   PointXYZ{(*icoords_)(idx[2], size_t(0)), (*icoords_)(idx[2], size_t(1)),
                                            (*icoords_)(idx[2], size_t(2))});

            // pick an epsilon based on a characteristic length (sqrt(area))
            // (this scales linearly so it better compares with linear weights u,v,w)
            const double edgeEpsilon = parametricEpsilon * std::sqrt(triag.area());
            ATLAS_ASSERT(edgeEpsilon >= 0);

            Intersect is = triag.intersects(ray, edgeEpsilon);

            if (is) {
                // weights are the linear Lagrange function evaluated at u,v (aka
                // barycentric coordinates)
                w[0] = 1. - is.u - is.v;
                w[1] = is.u;
                w[2] = is.v;

                if (on_triag_edge()) {
                    if (on_single_point()) {
                        triplets.emplace_back(ip, idx[single_point], w[single_point]);
                    }
                    else {
                        if ((*igidx_)(idx[edge.idx[1]]) < (*igidx_)(idx[edge.idx[0]])) {
                            edge.swap();
                        }
                        interpolate_edge(triag.p(edge.idx[0]), triag.p(edge.idx[1]));
                        for (size_t i = 0; i < 2; ++i) {
                            triplets.emplace_back(ip, idx[edge.idx[i]], w[edge.idx[i]]);
                        }
                    }
                }
                else {
                    for (size_t i = 0; i < 3; ++i) {
                        triplets.emplace_back(ip, idx[i], w[i]);
                    }
                }

                break;  // stop looking for elements
            }
        }
        else {
            /* quadrilateral */
            element::Quad3D quad(PointXYZ{(*icoords_)(idx[0], (size_t)0), (*icoords_)(idx[0], (size_t)1),
                                          (*icoords_)(idx[0], (size_t)2)},
                                 PointXYZ{(*icoords_)(idx[1], (size_t)0), (*icoords_)(idx[1], (size_t)1),
                                          (*icoords_)(idx[1], (size_t)2)},
                                 PointXYZ{(*icoords_)(idx[2], (size_t)0), (*icoords_)(idx[2], (size_t)1),
                                          (*icoords_)(idx[2], (size_t)2)},
                                 PointXYZ{(*icoords_)(idx[3], (size_t)0), (*icoords_)(idx[3], (size_t)1),
                                          (*icoords_)(idx[3], (size_t)2)});

            // pick an epsilon based on a characteristic length (sqrt(area))
            // (this scales linearly so it better compares with linear weights u,v,w)
            const double edgeEpsilon = parametricEpsilon * std::sqrt(quad.area());
            ATLAS_ASSERT(edgeEpsilon >= 0);

            Intersect is = quad.intersects(ray, edgeEpsilon);

            if (is) {
                // weights are the bilinear Lagrange function evaluated at u,v
                w[0] = (1. - is.u) * (1. - is.v);
                w[1] = is.u * (1. - is.v);
                w[2] = is.u * is.v;
                w[3] = (1. - is.u) * is.v;

                if (on_quad_edge()) {
                    if (on_single_point()) {
                        triplets.emplace_back(ip, idx[single_point], w[single_point]);
                    }
                    else {
                        if ((*igidx_)(idx[edge.idx[1]]) < (*igidx_)(idx[edge.idx[0]])) {
                            edge.swap();
                        }
                        interpolate_edge(quad.p(edge.idx[0]), quad.p(edge.idx[1]));
                        for (size_t i = 0; i < 2; ++i) {
                            triplets.emplace_back(ip, idx[edge.idx[i]], w[edge.idx[i]]);
                        }
                    }
                }
                else {
                    for (size_t i = 0; i < 4; ++i) {
                        triplets.emplace_back(ip, idx[i], w[i]);
                    }
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
