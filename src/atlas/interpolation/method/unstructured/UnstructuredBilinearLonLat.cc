/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iomanip>
#include <limits>

#include "UnstructuredBilinearLonLat.h"

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

MethodBuilder<UnstructuredBilinearLonLat> __builder_2("unstructured-bilinear-lonlat");

// epsilon used to scale edge tolerance when projecting ray to intesect element
static const double parametricEpsilon = 1e-15;

}  // namespace


void UnstructuredBilinearLonLat::do_setup(const Grid& source, const Grid& target, const Cache& cache) {
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

void UnstructuredBilinearLonLat::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    ATLAS_TRACE("atlas::interpolation::method::BilinearRemapping::do_setup()");

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

void UnstructuredBilinearLonLat::print(std::ostream& out) const {
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

void UnstructuredBilinearLonLat::setup(const FunctionSpace& source) {
    const functionspace::NodeColumns src = source;
    ATLAS_ASSERT(src);

    Mesh meshSource = src.mesh();


    auto trace_setup_source = atlas::Trace{Here(), "Setup source"};

    // 3D point coordinates
    Field source_xyz = mesh::actions::BuildXYZField("xyz")(meshSource);

    // generate barycenters of each triangle & insert them on a kd-tree
    util::Config config;
    config.set("name", "centre ");
    config.set("flatten_virtual_elements", false);
    Field cell_centres = mesh::actions::BuildCellCentres(config)(meshSource);

    std::unique_ptr<ElemIndex3> eTree(create_element_kdtree(meshSource, cell_centres));

    trace_setup_source.stop();

    ilonlat_.reset(new array::ArrayView<double, 2>(array::make_view<double, 2>(meshSource.nodes().lonlat())));
    olonlat_.reset(new array::ArrayView<double, 2>(array::make_view<double, 2>(target_lonlat_)));
    oxyz_.reset(new array::ArrayView<double, 2>(array::make_view<double, 2>(target_xyz_)));
    igidx_.reset(new array::ArrayView<gidx_t, 1>(array::make_view<gidx_t, 1>(src.nodes().global_index())));
    connectivity_              = &meshSource.cells().node_connectivity();
    const mesh::Nodes& i_nodes = meshSource.nodes();


    idx_t inp_npts = i_nodes.size();
    idx_t out_npts = olonlat_->shape(0);

    // return early if no output points on this partition reserve is called on
    // the triplets but also during the sparseMatrix constructor. This won't
    // work for empty matrices
    if (out_npts == 0) {
        return;
    }

    array::ArrayView<int, 1> out_ghosts = array::make_view<int, 1>(target_ghost_);

    idx_t Nelements = meshSource.cells().size();

    // weights -- one per vertex of element, triangles (3) or quads (4)

    Triplets weights_triplets;               // structure to fill-in sparse matrix
    weights_triplets.reserve(out_npts * 4);  // preallocate space as if all elements where quads

    // search nearest k cell centres

    const idx_t maxNbElemsToTry = std::max<idx_t>(8, idx_t(Nelements * max_fraction_elems_to_try_));

    std::vector<idx_t> failures;

    ATLAS_TRACE_SCOPE("Computing interpolation matrix") {
        eckit::ProgressTimer progress("Computing interpolation weights", out_npts, "point", double(5), Log::debug());
        for (idx_t ip = 0; ip < out_npts; ++ip, ++progress) {
            if (out_ghosts(ip)) {
                continue;
            }

            PointXYZ p{(*oxyz_)(ip, XX), (*oxyz_)(ip, YY), (*oxyz_)(ip, ZZ)};  // lookup point

            idx_t kpts   = 8;
            bool success = false;
            std::ostringstream failures_log;

            ElemIndex3::NodeList cs = eTree->kNearestNeighbours(p, kpts);
            Triplets triplets       = projectPointToElements(ip, cs, failures_log);

            if (triplets.size()) {
                std::copy(triplets.begin(), triplets.end(), std::back_inserter(weights_triplets));
                success = true;
            }

            if (!success) {
                failures.push_back(ip);
                if (not treat_failure_as_missing_value_) {
                    Log::debug() << "------------------------------------------------------"
                                    "---------------------\n";
                    Log::debug() << "Failed to project point (lon,lat)=" << (*olonlat_)(ip, LON) << " "
                                 << (*olonlat_)(ip, LAT) << '\n';
                    Log::debug() << failures_log.str();
                }
            }
        }
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
            for (auto i : failures) {
                const PointLonLat pll{(*olonlat_)(i, LON), (*olonlat_)(i, LAT)};  // lookup point
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

Method::Triplets UnstructuredBilinearLonLat::projectPointToElements(size_t ip, const ElemIndex3::NodeList& elems,
                                                                    std::ostream& /* failures_log */) const {
    ATLAS_ASSERT(elems.begin() != elems.end());

    const size_t inp_points = ilonlat_->shape(0);
    std::array<size_t, 4> idx;
    std::array<double, 4> w;
    std::array<double, 4> inv_dist_w;

    Triplets triplets;
    triplets.reserve(4);

    double o_lon{(*olonlat_)(ip, LON)};
    while (o_lon >= 360.0) {
        o_lon -= 360.0;
    }
    while (o_lon < 0.0) {
        o_lon += 360.0;
    }
    PointLonLat o_loc{o_lon, (*olonlat_)(ip, LAT)};  // lookup point

    auto inv_dist_weight_quad = [](element::Quad2D& q, const PointXY& loc, std::array<double, 4>& w) {
        double d[4];
        d[0] = util::Earth::distance({q.p(0).data()}, loc);
        d[1] = util::Earth::distance({q.p(1).data()}, loc);
        d[2] = util::Earth::distance({q.p(2).data()}, loc);
        d[3] = util::Earth::distance({q.p(3).data()}, loc);
        w[0] = d[1] * d[2] * d[3];
        w[1] = d[0] * d[2] * d[3];
        w[2] = d[1] * d[0] * d[3];
        w[3] = d[1] * d[0] * d[2];

        double suminv = 1. / (w[0] + w[1] + w[2] + w[3]);
        for (size_t i = 0; i < 4; ++i) {
            w[i] *= suminv;
        }
    };
    auto inv_dist_weight_triag = [](element::Triag2D& q, const PointXY& loc, std::array<double, 4>& w) {
        double d[3];
        d[0] = util::Earth::distance({q.p(0).data()}, loc);
        d[1] = util::Earth::distance({q.p(1).data()}, loc);
        d[2] = util::Earth::distance({q.p(2).data()}, loc);
        w[0] = d[1] * d[2];
        w[1] = d[0] * d[2];
        w[2] = d[1] * d[0];
        w[3] = 0.;

        double suminv = 1. / (w[0] + w[1] + w[2]);
        for (size_t i = 0; i < 3; ++i) {
            w[i] *= suminv;
        }
    };


    for (ElemIndex3::NodeList::const_iterator itc = elems.begin(); itc != elems.end(); ++itc) {
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
            element::Triag2D triag(PointLonLat{(*ilonlat_)(idx[0], LON), (*ilonlat_)(idx[0], LAT)},
                                   PointLonLat{(*ilonlat_)(idx[1], LON), (*ilonlat_)(idx[1], LAT)},
                                   PointLonLat{(*ilonlat_)(idx[2], LON), (*ilonlat_)(idx[2], LAT)});

            if (itc == elems.begin()) {
                inv_dist_weight_triag(triag, o_loc, inv_dist_w);
            }

            // pick an epsilon based on a characteristic length (sqrt(area))
            // (this scales linearly so it better compares with linear weights u,v,w)
            const double edgeEpsilon = parametricEpsilon * std::sqrt(triag.area());
            ATLAS_ASSERT(edgeEpsilon >= 0);

            Intersect is = triag.intersects(o_loc, edgeEpsilon);

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
            double lons[4];
            lons[0] = (*ilonlat_)(idx[0], LON);
            lons[1] = (*ilonlat_)(idx[1], LON);
            lons[2] = (*ilonlat_)(idx[2], LON);
            lons[3] = (*ilonlat_)(idx[3], LON);
            // adjust quadrilaterals to lie within [0,360]
            for (idx_t i = 0; i < 4; ++i) {
                while (lons[i] > 360.0) {
                    lons[i] -= 360.0;
                }
                while (lons[i] < 0.0) {
                    lons[i] += 360.0;
                }
            }
            // shift cells on the east-west periodic boundary from the east to the west
            // so that the quad surrounds a point with output longitude in [0,360)
            double minlon = std::numeric_limits<double>::max();
            for ( int i = 0; i < 4; i++ ) {
                minlon = std::min( minlon, lons[i] );
            }
            for ( int i = 0; i < 4; i++ ) {
                if ( (lons[i] - minlon) > 180 ) {
                    lons[i] -= 360;
                }
            }

            element::Quad2D quad(
                PointLonLat{lons[0], (*ilonlat_)(idx[0], LAT)}, PointLonLat{lons[1], (*ilonlat_)(idx[1], LAT)},
                PointLonLat{lons[2], (*ilonlat_)(idx[2], LAT)}, PointLonLat{lons[3], (*ilonlat_)(idx[3], LAT)});

            ATLAS_ASSERT( quad.validate() );

            if (itc == elems.begin()) {
                inv_dist_weight_quad(quad, o_loc, inv_dist_w);
            }

            // pick an epsilon based on a characteristic length (sqrt(area))
            // (this scales linearly so it better compares with linear weights u,v,w)
            const double edgeEpsilon = parametricEpsilon * std::sqrt(quad.area());
            ATLAS_ASSERT(edgeEpsilon >= 0);

            Intersect is = quad.localRemap(o_loc, edgeEpsilon);

            if (!is) {
                // repeat the calculation to catch points which are near the east boundary of the grid.
                is = quad.localRemap({o_lon - 360, (*olonlat_)(ip, LAT)}, edgeEpsilon);
            }

            if (is) {
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

    if (triplets.empty()) {
        // Crude inverse distance weighting to catch cells that are difficult
        // to identify and interpolate in a lon/lat projection, near the north
        // pole
        const idx_t elem_id = idx_t((*elems.begin()).value().payload());
        for (size_t i = 0; i < connectivity_->cols(elem_id); ++i) {
            idx[i] = (*connectivity_)(elem_id, i);
            triplets.emplace_back(ip, idx[i], inv_dist_w[i]);
        }
    }

    if (!triplets.empty()) {
        normalise(triplets);
    }
    return triplets;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
