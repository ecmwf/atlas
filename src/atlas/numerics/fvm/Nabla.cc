/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/config/Parametrisation.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/numerics/fvm/Nabla.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

// =======================================================

using Topology = atlas::mesh::Nodes::Topology;
using Range    = atlas::array::Range;

namespace atlas {
namespace numerics {
namespace fvm {

namespace {
static NablaBuilder<Nabla> __fvm_nabla("fvm");
}

Nabla::Nabla(const numerics::Method& method, const eckit::Parametrisation& p): atlas::numerics::NablaImpl(method, p) {
    fvm_ = dynamic_cast<const fvm::Method*>(&method);
    if (!fvm_) {
        throw_Exception("atlas::numerics::fvm::Nabla needs a atlas::numerics::fvm::Method", Here());
    }
    Log::debug() << "Nabla constructed for method " << fvm_->name() << " with "
                 << fvm_->node_columns().nb_nodes_global() << " nodes total" << std::endl;
    fvm_->attach();

    p.get("metric_approach", metric_approach_);
    // Experimental, only affects divergence and curl !!!
    // metric_approach = 0  ORIGINAL, DEFAULT
    //   metric_term cos(y) is multiplied with each wind component
    //   --> cell-interface:  avg = 0.5*( cos(y1)*u1 + cos(y2)*u2 )
    // metric_approach = 1:
    //   metric_term cos(0.5*(y1+y2)) is used at cell interface
    //   --> cell-interface: avg = 0.5*(u1+u2)*cos(0.5*(y1+y2))
    // Results seem to indicate that approach=0 is overall better, although approach=1
    // seems to handle pole slightly better (error factor 2 to 4 times lower)

    setup();
}

Nabla::~Nabla() = default;

void Nabla::setup() {
    const mesh::Edges& edges = fvm_->mesh().edges();

    const idx_t nedges = fvm_->edge_columns().nb_edges();

    const auto edge_flags = array::make_view<int, 1>(edges.flags());
    auto is_pole_edge     = [&](idx_t e) { return Topology::check(edge_flags(e), Topology::POLE); };

    // Filter pole_edges out of all edges
    std::vector<idx_t> tmp(nedges);
    idx_t c(0);
    for (idx_t jedge = 0; jedge < nedges; ++jedge) {
        if (is_pole_edge(jedge)) {
            tmp[c++] = jedge;
        }
    }
    pole_edges_.clear();
    pole_edges_.reserve(c);
    for (idx_t jedge = 0; jedge < c; ++jedge) {
        pole_edges_.push_back(tmp[jedge]);
    }
}

void Nabla::gradient(const Field& field, Field& grad_field) const {
    if (field.variables() > 1) {
        return gradient_of_vector(field, grad_field);
    }
    else {
        return gradient_of_scalar(field, grad_field);
    }
}

void Nabla::gradient_of_scalar(const Field& scalar_field, Field& grad_field) const {
    Log::debug() << "Compute gradient of scalar field " << scalar_field.name() << " with fvm method" << std::endl;

    auto dispatch = [&](auto value) {
        using Value          = std::decay_t<decltype(value)>;

        const Value radius  = fvm_->radius();
        const Value deg2rad = M_PI / 180.;

        const mesh::Edges& edges = fvm_->mesh().edges();
        const mesh::Nodes& nodes = fvm_->mesh().nodes();

        const idx_t nnodes = fvm_->node_columns().nb_nodes();
        const idx_t nedges = fvm_->edge_columns().nb_edges();

        const auto scalar = scalar_field.levels()
                                ? array::make_view<Value, 2>(scalar_field).slice(Range::all(), Range::all())
                                : array::make_view<Value, 1>(scalar_field).slice(Range::all(), Range::dummy());
        auto grad         = grad_field.levels()
                                ? array::make_view<Value, 3>(grad_field).slice(Range::all(), Range::all(), Range::all())
                                : array::make_view<Value, 2>(grad_field).slice(Range::all(), Range::dummy(), Range::all());

        const idx_t nlev = scalar.shape(1);
        if (grad.shape(1) != nlev) {
            throw_AssertionFailed("gradient field should have same number of levels", Here());
        }

        const auto lonlat_deg     = array::make_view<double, 2>(nodes.lonlat());
        const auto dual_volumes   = array::make_view<double, 1>(nodes.field("dual_volumes"));
        const auto dual_normals   = array::make_view<double, 2>(edges.field("dual_normals"));
        const auto node2edge_sign = array::make_view<double, 2>(nodes.field("node2edge_sign"));

        const mesh::Connectivity& node2edge           = nodes.edge_connectivity();
        const mesh::MultiBlockConnectivity& edge2node = edges.node_connectivity();

        array::ArrayT<Value> avgS_arr(nedges, nlev, 2ul);
        auto avgS = array::make_view<Value, 3>(avgS_arr);

        const Value scale = deg2rad * deg2rad * radius;

        atlas_omp_parallel {
            atlas_omp_for(idx_t jedge = 0; jedge < nedges; ++jedge) {
                idx_t ip1 = edge2node(jedge, 0);
                idx_t ip2 = edge2node(jedge, 1);
                Value S[2] = { static_cast<Value>(dual_normals(jedge, LON)), static_cast<Value>(dual_normals(jedge, LAT))};

                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    Value avg              = (scalar(ip1, jlev) + scalar(ip2, jlev)) * Value{0.5};
                    avgS(jedge, jlev, LON) = S[LON] * deg2rad * avg;
                    avgS(jedge, jlev, LAT) = S[LAT] * deg2rad * avg;
                }
            }

            atlas_omp_for(idx_t jnode = 0; jnode < nnodes; ++jnode) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    grad(jnode, jlev, LON) = 0.;
                    grad(jnode, jlev, LAT) = 0.;
                }
                for (idx_t jedge = 0; jedge < node2edge.cols(jnode); ++jedge) {
                    const idx_t iedge = node2edge(jnode, jedge);
                    if (iedge < nedges) {
                        const Value add = node2edge_sign(jnode, jedge);
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            grad(jnode, jlev, LON) += add * avgS(iedge, jlev, LON);
                            grad(jnode, jlev, LAT) += add * avgS(iedge, jlev, LAT);
                        }
                    }
                }

                const Value y        = lonlat_deg(jnode, LAT) * deg2rad;
                const Value metric_y = Value{1.} / (static_cast<Value>(dual_volumes(jnode)) * scale);
                const Value metric_x = metric_y / std::cos(y);
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    grad(jnode, jlev, LON) *= metric_x;
                    grad(jnode, jlev, LAT) *= metric_y;
                }
            }
        }
    };
    ATLAS_ASSERT( scalar_field.datatype() == grad_field.datatype() );
    switch (scalar_field.datatype().kind()) {
        case (DataType::KIND_REAL32): {
            dispatch(float{});
            break;
        }
        case (DataType::KIND_REAL64): {
            dispatch(double{});
            break;
        }
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

// ================================================================================

void Nabla::gradient_of_vector(const Field& vector_field, Field& grad_field) const {
    Log::debug() << "Compute gradient of vector field " << vector_field.name() << " with fvm method" << std::endl;

    auto dispatch = [&](auto value) {
        using Value          = std::decay_t<decltype(value)>;

        const Value radius  = fvm_->radius();
        const Value deg2rad = M_PI / 180.;

        const mesh::Edges& edges = fvm_->mesh().edges();
        const mesh::Nodes& nodes = fvm_->mesh().nodes();

        const idx_t nnodes = fvm_->node_columns().nb_nodes();
        const idx_t nedges = fvm_->edge_columns().nb_edges();

        const auto vector =
            vector_field.levels()
                ? array::make_view<Value, 3>(vector_field).slice(Range::all(), Range::all(), Range::all())
                : array::make_view<Value, 2>(vector_field).slice(Range::all(), Range::dummy(), Range::all());
        auto grad = grad_field.levels()
                        ? array::make_view<Value, 3>(grad_field).slice(Range::all(), Range::all(), Range::all())
                        : array::make_view<Value, 2>(grad_field).slice(Range::all(), Range::dummy(), Range::all());

        const idx_t nlev = vector.shape(1);
        if (grad.shape(1) != nlev) {
            throw_AssertionFailed("gradient field should have same number of levels", Here());
        }

        const auto lonlat_deg     = array::make_view<double, 2>(nodes.lonlat());
        const auto dual_volumes   = array::make_view<double, 1>(nodes.field("dual_volumes"));
        const auto dual_normals   = array::make_view<double, 2>(edges.field("dual_normals"));
        const auto node2edge_sign = array::make_view<double, 2>(nodes.field("node2edge_sign"));
        const auto edge_flags     = array::make_view<int, 1>(edges.flags());
        auto is_pole_edge         = [&](idx_t e) { return Topology::check(edge_flags(e), Topology::POLE); };

        const mesh::Connectivity& node2edge           = nodes.edge_connectivity();
        const mesh::MultiBlockConnectivity& edge2node = edges.node_connectivity();

        array::ArrayT<Value> avgS_arr(nedges, nlev, 4ul);
        array::ArrayView<Value, 3> avgS = array::make_view<Value, 3>(avgS_arr);

        const Value scale = deg2rad * deg2rad * radius;

        enum
        {
            LONdLON = 0,
            LONdLAT = 1,
            LATdLON = 2,
            LATdLAT = 3
        };

        atlas_omp_parallel {
            atlas_omp_for(idx_t jedge = 0; jedge < nedges; ++jedge) {
                idx_t ip1  = edge2node(jedge, 0);
                idx_t ip2  = edge2node(jedge, 1);
                Value pbc = 1. - 2. * is_pole_edge(jedge);

                Value S[2] = {static_cast<Value>(dual_normals(jedge,LON)), static_cast<Value>(dual_normals(jedge,LAT))};
                Value avg[2];

                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    avg[LON] = (vector(ip1, jlev, LON) + pbc * vector(ip2, jlev, LON)) * Value{0.5};
                    avg[LAT] = (vector(ip1, jlev, LAT) + pbc * vector(ip2, jlev, LAT)) * Value{0.5};
                    avgS(jedge, jlev, LONdLON) = S[LON] * deg2rad * avg[LON];
                    // above = 0 at pole because of dual_normals
                    avgS(jedge, jlev, LONdLAT) = S[LAT] * deg2rad * avg[LON];
                    avgS(jedge, jlev, LATdLON) = S[LON] * deg2rad * avg[LAT];
                    // above = 0 at pole because of dual_normals
                    avgS(jedge, jlev, LATdLAT) = S[LAT] * deg2rad * avg[LAT];
                }
            }

            atlas_omp_for(idx_t jnode = 0; jnode < nnodes; ++jnode) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    grad(jnode, jlev, LONdLON) = 0.;
                    grad(jnode, jlev, LONdLAT) = 0.;
                    grad(jnode, jlev, LATdLON) = 0.;
                    grad(jnode, jlev, LATdLAT) = 0.;
                }
                for (idx_t jedge = 0; jedge < node2edge.cols(jnode); ++jedge) {
                    const idx_t iedge = node2edge(jnode, jedge);
                    if (iedge < nedges) {
                        Value add = node2edge_sign(jnode, jedge);
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            grad(jnode, jlev, LONdLON) += add * avgS(iedge, jlev, LONdLON);
                            grad(jnode, jlev, LONdLAT) += add * avgS(iedge, jlev, LONdLAT);
                            grad(jnode, jlev, LATdLON) += add * avgS(iedge, jlev, LATdLON);
                            grad(jnode, jlev, LATdLAT) += add * avgS(iedge, jlev, LATdLAT);
                        }
                    }
                }
                const Value y        = lonlat_deg(jnode, LAT) * deg2rad;
                const Value metric_y = Value{1.} / (static_cast<Value>(dual_volumes(jnode)) * scale);
                const Value metric_x = metric_y / std::cos(y);
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    grad(jnode, jlev, LONdLON) *= metric_x;
                    grad(jnode, jlev, LATdLON) *= metric_x;
                    grad(jnode, jlev, LONdLAT) *= metric_y;
                    grad(jnode, jlev, LATdLAT) *= metric_y;
                }
            }
        }
        // Fix wrong node2edge_sign for vector quantities
        for (size_t jedge = 0; jedge < pole_edges_.size(); ++jedge) {
            const idx_t iedge     = pole_edges_[jedge];
            const idx_t jnode     = edge2node(iedge, 1);
            const Value metric_y = Value{1.} / (static_cast<Value>(dual_volumes(jnode)) * scale);
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                grad(jnode, jlev, LONdLAT) -= 2. * avgS(iedge, jlev, LONdLAT) * metric_y;
                grad(jnode, jlev, LATdLAT) -= 2. * avgS(iedge, jlev, LATdLAT) * metric_y;
            }
        }
    };
    ATLAS_ASSERT( vector_field.datatype() == grad_field.datatype() );
    switch (vector_field.datatype().kind()) {
        case (DataType::KIND_REAL32): {
            dispatch(float{});
            break;
        }
        case (DataType::KIND_REAL64): {
            dispatch(double{});
            break;
        }
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

// ================================================================================

void Nabla::divergence(const Field& vector_field, Field& div_field) const {

    auto dispatch = [&](auto value) {
        using Value          = std::decay_t<decltype(value)>;

        const Value radius  = fvm_->radius();
        const Value deg2rad = M_PI / 180.;

        const mesh::Edges& edges = fvm_->mesh().edges();
        const mesh::Nodes& nodes = fvm_->mesh().nodes();

        const idx_t nnodes = fvm_->node_columns().nb_nodes();
        const idx_t nedges = fvm_->edge_columns().nb_edges();

        const auto vector =
            vector_field.levels()
                ? array::make_view<Value, 3>(vector_field).slice(Range::all(), Range::all(), Range::all())
                : array::make_view<Value, 2>(vector_field).slice(Range::all(), Range::dummy(), Range::all());
        auto div = div_field.levels() ? array::make_view<Value, 2>(div_field).slice(Range::all(), Range::all())
                                      : array::make_view<Value, 1>(div_field).slice(Range::all(), Range::dummy());

        const idx_t nlev = vector.shape(1);
        if (div.shape(1) != nlev) {
            throw_AssertionFailed("div_field should have same number of levels", Here());
        }

        const auto lonlat_deg     = array::make_view<double, 2>(nodes.lonlat());
        const auto dual_volumes   = array::make_view<double, 1>(nodes.field("dual_volumes"));
        const auto dual_normals   = array::make_view<double, 2>(edges.field("dual_normals"));
        const auto node2edge_sign = array::make_view<double, 2>(nodes.field("node2edge_sign"));
        const auto edge_flags     = array::make_view<int, 1>(edges.flags());
        auto is_pole_edge         = [&](idx_t e) { return Topology::check(edge_flags(e), Topology::POLE); };

        const mesh::Connectivity& node2edge           = nodes.edge_connectivity();
        const mesh::MultiBlockConnectivity& edge2node = edges.node_connectivity();

        array::ArrayT<Value> avgS_arr(nedges, nlev, 2ul);
        array::ArrayView<Value, 3> avgS = array::make_view<Value, 3>(avgS_arr);

        const Value scale = deg2rad * deg2rad * radius;

        enum
        {
            LONdLON = 0,
            LATdLAT = 1
        };
        atlas_omp_parallel {
            atlas_omp_for(idx_t jedge = 0; jedge < nedges; ++jedge) {
                Value pbc = 1 - is_pole_edge(jedge);

                idx_t ip1 = edge2node(jedge, 0);
                idx_t ip2 = edge2node(jedge, 1);
                Value y1 = lonlat_deg(ip1, LAT) * deg2rad;
                Value y2 = lonlat_deg(ip2, LAT) * deg2rad;

                Value cosy1, cosy2;
                if (metric_approach_ == 0) {
                    cosy1 = std::cos(y1) * pbc;
                    cosy2 = std::cos(y2) * pbc;
                }
                else {
                    cosy1 = cosy2 = std::cos(0.5 * (y1 + y2)) * pbc;
                }

                Value S[2] = {static_cast<Value>(dual_normals(jedge, LON)) * deg2rad, static_cast<Value>(dual_normals(jedge, LAT)) * deg2rad};
                Value avg[2];

                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    Value u1 = vector(ip1, jlev, LON);
                    Value u2 = vector(ip2, jlev, LON);
                    Value v1 = vector(ip1, jlev, LAT) * cosy1;
                    Value v2 = vector(ip2, jlev, LAT) * cosy2;

                    avg[LON] = (u1 + u2) * 0.5;
                    avg[LAT] = (v1 + v2) * 0.5;

                    avgS(jedge, jlev, LONdLON) = avg[LON] * S[LON];
                    avgS(jedge, jlev, LATdLAT) = avg[LAT] * S[LAT];
                }
            }

            atlas_omp_for(idx_t jnode = 0; jnode < nnodes; ++jnode) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    div(jnode, jlev) = 0.;
                }
                for (idx_t jedge = 0; jedge < node2edge.cols(jnode); ++jedge) {
                    idx_t iedge = node2edge(jnode, jedge);
                    if (iedge < nedges) {
                        Value add = node2edge_sign(jnode, jedge);
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            div(jnode, jlev) += add * (avgS(iedge, jlev, LONdLON) + avgS(iedge, jlev, LATdLAT));
                        }
                    }
                }
                const Value y = lonlat_deg(jnode, LAT) * deg2rad;
                Value metric  = Value{1.} / (static_cast<Value>(dual_volumes(jnode)) * scale * std::cos(y));
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    div(jnode, jlev) *= metric;
                }
            }
        }

    };
    ATLAS_ASSERT( vector_field.datatype() == div_field.datatype() );
    switch (vector_field.datatype().kind()) {
        case (DataType::KIND_REAL32): {
            dispatch(float{});
            break;
        }
        case (DataType::KIND_REAL64): {
            dispatch(double{});
            break;
        }
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

void Nabla::curl(const Field& vector_field, Field& curl_field) const {
    auto dispatch = [&](auto value) {
        using Value          = std::decay_t<decltype(value)>;


        const Value radius  = fvm_->radius();
        const Value deg2rad = M_PI / 180.;

        const mesh::Edges& edges = fvm_->mesh().edges();
        const mesh::Nodes& nodes = fvm_->mesh().nodes();

        const idx_t nnodes = fvm_->node_columns().nb_nodes();
        const idx_t nedges = fvm_->edge_columns().nb_edges();

        const auto vector =
            vector_field.levels()
                ? array::make_view<Value, 3>(vector_field).slice(Range::all(), Range::all(), Range::all())
                : array::make_view<Value, 2>(vector_field).slice(Range::all(), Range::dummy(), Range::all());
        auto curl = curl_field.levels() ? array::make_view<Value, 2>(curl_field).slice(Range::all(), Range::all())
                                        : array::make_view<Value, 1>(curl_field).slice(Range::all(), Range::dummy());

        const idx_t nlev = vector.shape(1);
        if (curl.shape(1) != nlev) {
            throw_AssertionFailed("curl field should have same number of levels", Here());
        }

        const auto lonlat_deg     = array::make_view<double, 2>(nodes.lonlat());
        const auto dual_volumes   = array::make_view<double, 1>(nodes.field("dual_volumes"));
        const auto dual_normals   = array::make_view<double, 2>(edges.field("dual_normals"));
        const auto node2edge_sign = array::make_view<double, 2>(nodes.field("node2edge_sign"));
        const auto edge_flags     = array::make_view<int, 1>(edges.flags());
        auto is_pole_edge         = [&](idx_t e) { return Topology::check(edge_flags(e), Topology::POLE); };


        const mesh::Connectivity& node2edge           = nodes.edge_connectivity();
        const mesh::MultiBlockConnectivity& edge2node = edges.node_connectivity();

        array::ArrayT<Value> avgS_arr(nedges, nlev, 2ul);
        array::ArrayView<Value, 3> avgS = array::make_view<Value, 3>(avgS_arr);

        const Value scale = deg2rad * deg2rad * radius;

        enum
        {
            LONdLAT = 0,
            LATdLON = 1
        };


        atlas_omp_parallel {
            atlas_omp_for(idx_t jedge = 0; jedge < nedges; ++jedge) {
                idx_t ip1 = edge2node(jedge, 0);
                idx_t ip2 = edge2node(jedge, 1);
                Value y1 = lonlat_deg(ip1, LAT) * deg2rad;
                Value y2 = lonlat_deg(ip2, LAT) * deg2rad;

                Value pbc = 1 - is_pole_edge(jedge);
                Value cosy1;
                Value cosy2;
                if (metric_approach_ == 0) {
                    cosy1 = std::cos(y1) * pbc;
                    cosy2 = std::cos(y2) * pbc;
                }
                else {
                    cosy1 = cosy2 = std::cos(Value{0.5} * (y1 + y2)) * pbc;
                }

                Value S[2] = {static_cast<Value>(dual_normals(jedge, LON)) * deg2rad, static_cast<Value>(dual_normals(jedge, LAT)) * deg2rad};
                Value avg[2];

                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    Value u1 = vector(ip1, jlev, LON) * cosy1;
                    Value u2 = vector(ip2, jlev, LON) * cosy2;
                    Value v1 = vector(ip1, jlev, LAT);
                    Value v2 = vector(ip2, jlev, LAT);

                    avg[LON] = (u1 + u2) * Value{0.5};
                    avg[LAT] = (v1 + v2) * Value{0.5};

                    avgS(jedge, jlev, LONdLAT) = avg[LON] * S[LAT];
                    avgS(jedge, jlev, LATdLON) = avg[LAT] * S[LON];
                }
            }

            atlas_omp_for(idx_t jnode = 0; jnode < nnodes; ++jnode) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    curl(jnode, jlev) = 0.;
                }
                for (idx_t jedge = 0; jedge < node2edge.cols(jnode); ++jedge) {
                    idx_t iedge = node2edge(jnode, jedge);
                    if (iedge < nedges) {
                        double add = node2edge_sign(jnode, jedge);
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            curl(jnode, jlev) += add * (avgS(iedge, jlev, LATdLON) - avgS(iedge, jlev, LONdLAT));
                        }
                    }
                }
                Value y      = static_cast<Value>(lonlat_deg(jnode, LAT)) * deg2rad;
                Value metric = Value{1.} / (static_cast<Value>(dual_volumes(jnode)) * scale * std::cos(y));
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    curl(jnode, jlev) *= metric;
                }
            }
        }
    };
    ATLAS_ASSERT( vector_field.datatype() == curl_field.datatype() );
    switch (vector_field.datatype().kind()) {
        case (DataType::KIND_REAL32): {
            dispatch(float{});
            break;
        }
        case (DataType::KIND_REAL64): {
            dispatch(double{});
            break;
        }
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

void Nabla::laplacian(const Field& scalar, Field& lapl) const {
    Field grad(fvm_->node_columns().createField(option::name("grad") | option::levels(scalar.levels()) |
                                                option::variables(2) | option::datatype(scalar.datatype())));
    gradient(scalar, grad);
    if (fvm_->node_columns().halo().size() < 2) {
        fvm_->node_columns().haloExchange(grad);
    }
    divergence(grad, lapl);
}

const FunctionSpace& Nabla::functionspace() const {
    return fvm_->node_columns();
}

}  // namespace fvm
}  // namespace numerics
}  // namespace atlas
