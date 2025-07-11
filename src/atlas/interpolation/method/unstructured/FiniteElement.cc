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

void FiniteElement::do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache& cache) { 
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

void FiniteElement::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
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
            target_xyz_    = mesh::actions::BuildXYZField("xyz")(meshTarget);
        }
        else {
            target_xyz_    = create_xyz(target_lonlat_);
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

void FiniteElement::setup(const FunctionSpace& source) {
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
            std::copy(triplets.begin(), triplets.end(), weights_triplets.begin()+4*n);
            return true;
        }
        return false;
    };

    double search_radius = 0.;
    if (meshSource.metadata().has("cell_maximum_diagonal_on_unit_sphere")) {
        search_radius = Geometry("Earth").radius() * meshSource.metadata().getDouble("cell_maximum_diagonal_on_unit_sphere");
        ASSERT(search_radius > 0.);
        Log::debug() << "k-d tree: search radius = " << search_radius/1000. << " km" << std::endl;
    }
    auto find_element_candidates_in_search_radius = [&eTree,&search_radius](const PointXYZ& p) {
        return eTree->findInSphere(p, search_radius);
    };
    auto find_k_nearest_element_candidates = [&eTree](const PointXYZ& p, size_t k) {
        return eTree->kNearestNeighbours(p, k);
    };
    auto try_interpolate_with_element_candidates = [&insert_triplets, this](idx_t n, const ElemIndex3::NodeList& element_candidates) -> bool {
        if (element_candidates.empty()) {
            return false;
        }
        return insert_triplets(n, projectPointToElements(n, element_candidates));
    };


    // search nearest k cell centres
    const idx_t maxNbElemsToTry = std::max<idx_t>(8, idx_t(Nelements * max_fraction_elems_to_try_));
    size_t diagnosed_max_neighbours         = 0;
    bool allowed_to_diagnose_max_neighbours = Log::debug() && atlas_omp_get_max_threads() > 1;

    std::vector<size_t> failures;

    ATLAS_TRACE_SCOPE("Computing interpolation matrix") {
        std::unique_ptr<eckit::ProgressTimer> progress;
        if (atlas_omp_get_max_threads() == 1) {
            progress.reset(new eckit::ProgressTimer{"Computing interpolation weights", static_cast<size_t>(out_npts), "point", double(5), Log::debug()});
        }
        atlas_omp_parallel_for (idx_t ip = 0; ip < out_npts; ++ip) {
            if (out_ghosts(ip)) {
                continue;
            }

            bool success = false;
            if (search_radius != 0.) {
                auto p = target_point(ip);
                success = try_interpolate_with_element_candidates(ip, find_element_candidates_in_search_radius(p));
            }
            else {
                size_t k = 1;
                auto p = target_point(ip);
                while (!success && k <= maxNbElemsToTry) {
                    if (allowed_to_diagnose_max_neighbours) { // avoid race condition
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
        Log::debug() << "Maximum neighbours searched was " << eckit::Plural(diagnosed_max_neighbours, "element") << std::endl;
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

struct ElementEdge {
    std::array<idx_t, 2> idx;
    void swap() {
        idx_t tmp = idx[0];
        idx[0]    = idx[1];
        idx[1]    = tmp;
    }
};

Method::Triplets FiniteElement::projectPointToElements(size_t ip, const ElemIndex3::NodeList& elems) const {
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

        const idx_t nb_cols = [&]() {
            int nb_cols = connectivity_->cols(elem_id);
            if (nb_cols == 5) {
                // Check if pentagon degenerates to quad. Otherwise abort.
                // For now only check if first and last point coincide.
                auto i1    = (*connectivity_)(elem_id, 0);
                auto iN    = (*connectivity_)(elem_id, nb_cols - 1);
                auto first = PointXYZ{(*icoords_)(i1, XX), (*icoords_)(i1, YY), (*icoords_)(i1, ZZ)};
                auto last  = PointXYZ{(*icoords_)(iN, XX), (*icoords_)(iN, YY), (*icoords_)(iN, ZZ)};
                if (first == last) {
                    return 4;
                }
            }
            return nb_cols;
        }();

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
            element::Triag3D triag(PointXYZ{(*icoords_)(idx[0], XX), (*icoords_)(idx[0], YY), (*icoords_)(idx[0], ZZ)},
                                   PointXYZ{(*icoords_)(idx[1], XX), (*icoords_)(idx[1], YY), (*icoords_)(idx[1], ZZ)},
                                   PointXYZ{(*icoords_)(idx[2], XX), (*icoords_)(idx[2], YY), (*icoords_)(idx[2], ZZ)});

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
            element::Quad3D quad(PointXYZ{(*icoords_)(idx[0], XX), (*icoords_)(idx[0], YY), (*icoords_)(idx[0], ZZ)},
                                 PointXYZ{(*icoords_)(idx[1], XX), (*icoords_)(idx[1], YY), (*icoords_)(idx[1], ZZ)},
                                 PointXYZ{(*icoords_)(idx[2], XX), (*icoords_)(idx[2], YY), (*icoords_)(idx[2], ZZ)},
                                 PointXYZ{(*icoords_)(idx[3], XX), (*icoords_)(idx[3], YY), (*icoords_)(idx[3], ZZ)});

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
