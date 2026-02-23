/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <iomanip>
#include <unordered_set>
#include <vector>

#include "ConservativeSphericalPolygonInterpolationLimiter.h"

#include "atlas/grid.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/library/FloatingPointExceptions.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildNode2CellConnectivity.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Topology.h"
#include "atlas/util/detail/filesystem.h"

#include "eckit/log/Bytes.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/mpi/Comm.h"

#define PRINT_BAD_POLYGONS 0

namespace atlas {
namespace interpolation {
namespace method {


using runtime::trace::StopWatch;
using Polygon = util::ConvexSphericalPolygon;
using PolygonArray = std::vector<util::ConvexSphericalPolygon>;


ConservativeSphericalPolygonInterpolationLimiter::
ConservativeSphericalPolygonInterpolationLimiter(const ConservativeSphericalPolygonInterpolation& interpolation):
    interpolation_(interpolation), src_cell_data_(interpolation.src_cell_data_), tgt_cell_data_(interpolation.tgt_cell_data_),
    limiter_(interpolation.limiter_), order_(interpolation.order_),
    src_fs_(interpolation.source()), tgt_fs_(interpolation.target()), data_(interpolation.data_) {
    // sharable_data_ = std::make_shared<Data>();
    // cache_         = Cache(sharable_data_);
    // data_          = sharable_data_.get();
}


double ConservativeSphericalPolygonInterpolationLimiter::limit(const Field& src_field, Field& tgt_field) {
    double mass_change = 0.;
    if (order_ == 2 && (limiter_ != "none")) {
        // set this environment variable to replace target_field with values showing
        //   0: target_field as it is
        //   1: values corrections from the limiter
        //   2: tcells marked for fixing
        unsigned int limiter_override_tgt = 0;
        const char* ATLAS_INTERPOLATION_LIMITER = ::getenv("ATLAS_INTERPOLATION_LIMITER");
        if (ATLAS_INTERPOLATION_LIMITER != nullptr) {
                limiter_override_tgt = std::atof(ATLAS_INTERPOLATION_LIMITER);
        }
        if (limiter_override_tgt > 2) {
            Log::error() << "ATLAS_INTERPOLATION_LIMITER can be:\n";
            Log::error() << "\t0 (default)\n";
            Log::error() << "\t1 (limiter correction field)\n";
            Log::error() << "\t2 (violation & collateral target cells)" << std::endl;
            ATLAS_ASSERT(false);
        }
        Field tgt_lim_field = tgt_fs_.createField<double>();
        auto tgt_lim_vals   = array::make_view<double, 1>(tgt_lim_field);
        for (idx_t tcell = 0; tcell < tgt_lim_vals.size(); ++tcell) {
            tgt_lim_vals(tcell) = 0.;
        }

        struct SrcActed {
            Indices tcells_done;
        };
        std::vector<SrcActed> src_acted_tgt;
        const auto src_vals = array::make_view<double, 1>(src_field);
        auto tgt_vals       = array::make_view<double, 1>(tgt_field);
        const auto& tgt_iparam = data_->tgt_iparam_;
        const auto& src_iparam = data_->src_iparam_;
        const auto& src_points = data_->src_.points;
        const auto& tgt_areas = data_->tgt_.areas;
        std::vector<PointXYZ> src_grads;

        if (tgt_cell_data_ && src_cell_data_) {
            ATLAS_TRACE_SCOPE("Compute source gradients") {
                src_grads.resize(src_vals.size());
                for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
                    src_grads[scell] = interpolation_.src_gradient_celldata(scell, src_vals);
                }
            }
            src_acted_tgt.resize(src_vals.size());
            auto& mpi_comm = mpi::comm();
            std::set<Indices> send_marked_scells_set;
            ConservativeSphericalPolygonInterpolation::Workspace_get_cell_neighbours w_cell;

            for (idx_t tcsp = 0; tcsp < data_->tgt_.csp_size; ++tcsp) {
                const auto& iparam = tgt_iparam[tcsp];
                double smin = std::numeric_limits<double>::max();
                double smax = -smin;
                for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                    idx_t scsp_id = iparam.csp_ids[i_scsp];
                    idx_t scell   = interpolation_.csp_to_cell(scsp_id, data_->src_);
                    smax = std::max(smax, src_vals(scell));
                    smin = std::min(smin, src_vals(scell));
                    const auto src_neighbours = interpolation_.get_cell_neighbours(interpolation_.src_mesh_, scell, w_cell);
                    for (auto nb : src_neighbours) {
                        smax = std::max(smax, src_vals(nb));
                        smin = std::min(smin, src_vals(nb));
                    }
                }
                idx_t tcell = interpolation_.csp_to_cell(tcsp, data_->tgt_);
                bool undershoot = (tgt_vals(tcell) < smin);
                bool overshoot = (tgt_vals(tcell) > smax);
                if (undershoot || overshoot) {
                    if (limiter_ == "zeroslope") {
                        send_marked_scells_set.insert(iparam.csp_ids);

                        if (limiter_override_tgt == 2) {
                            tgt_lim_vals(tcell) = 1.;
                        }
                        for (idx_t i_scsp = 0; i_scsp < iparam.csp_ids.size(); ++i_scsp) {
                            idx_t scsp_id = iparam.csp_ids[i_scsp];
                            idx_t scell   = interpolation_.csp_to_cell(scsp_id, data_->src_);
                            const PointXYZ& src_barycentre = src_points[scell];
                            PointXYZ scell_grad  = src_grads[scell];
                            scell_grad           = scell_grad - PointXYZ::mul(src_barycentre, PointXYZ::dot(scell_grad, src_barycentre));
                            auto& siparam = src_iparam[scsp_id];
                            for (idx_t i_tcsp_collateral = 0; i_tcsp_collateral < siparam.csp_ids.size(); ++i_tcsp_collateral) {
                                auto tcsp_collateral = siparam.csp_ids[i_tcsp_collateral];
                                const auto& iparam_collateral = tgt_iparam[tcsp_collateral];
                                auto tcell_collateral = interpolation_.csp_to_cell(tcsp_collateral, data_->tgt_);
                                // find the index of scell entry in iparam_collateral.csp_ids
                                auto scell_it = std::find(iparam_collateral.csp_ids.begin(), iparam_collateral.csp_ids.end(), scell);
                                idx_t scell_idx = scell_it - iparam_collateral.csp_ids.begin();
                                ATLAS_ASSERT(iparam_collateral.csp_ids[scell_idx] == scell);
                                double tgt_lim_val = iparam_collateral.weights[scell_idx] * PointXYZ::dot(scell_grad, iparam_collateral.centroids[scell_idx] - src_barycentre);
                                if (tgt_areas[tcell_collateral] > 0.) {
                                    tgt_lim_val /= tgt_areas[tcell_collateral];
                                }
                                SrcActed& it = src_acted_tgt[scell];
                                if (std::find(it.tcells_done.begin(), it.tcells_done.end(), tcell_collateral) == it.tcells_done.end()) {
                                    it.tcells_done.push_back(tcell_collateral);
                                    if (limiter_override_tgt == 2 && tgt_lim_vals(tcell_collateral) < 0.5) {
                                        tgt_lim_vals(tcell_collateral) = -1;
                                    }
                                    else {
                                        tgt_lim_vals(tcell_collateral) -= tgt_lim_val;
                                    }
                                }
                            }
                        }
                    }
                    else if (limiter_ == "clip") {
                        if (overshoot) {
                            tgt_lim_vals(tcell) = smax - tgt_vals(tcell);
                        }
                        else if (undershoot) {
                            tgt_lim_vals(tcell) = smin - tgt_vals(tcell);
                        }
                    }
                }
            }

            std::vector<Indices> send_marked_scells(send_marked_scells_set.begin(), send_marked_scells_set.end());
            eckit::mpi::Buffer<Indices> recv_marked_scells_buf(mpi_comm.size());
            // mpi_comm.allGatherv(send_marked_scells.begin(), send_marked_scells.end(), recv_marked_scells_buf);

            if (! limiter_override_tgt) {
                if (interpolation_.matrix_free_) {
                    for (idx_t tcell = 0 ; tcell < tgt_vals.size(); ++tcell) {
                        mass_change += tgt_lim_vals(tcell)  * tgt_areas[tcell];
                        tgt_vals(tcell) += tgt_lim_vals(tcell);
                    }
                }
                else {
                    for (idx_t tcell = 0 ; tcell < tgt_vals.size(); ++tcell) {
                        mass_change += tgt_lim_vals(tcell) * tgt_areas[tcell];
                        tgt_vals(tcell) -= tgt_lim_vals(tcell);
                    }
                }
            }
            else {
                for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                    mass_change += (tgt_lim_vals(tcell) - tgt_vals(tcell)) * tgt_areas[tcell];
                    tgt_vals(tcell) = tgt_lim_vals(tcell);
                }
            }
        }
        else {
            Log::info() << "Limiter supports only CellColumns data." << std::endl;
            ATLAS_NOTIMPLEMENTED;
        }
    }
    return mass_change;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas