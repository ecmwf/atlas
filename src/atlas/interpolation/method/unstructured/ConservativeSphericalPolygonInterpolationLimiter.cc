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


namespace atlas {
namespace interpolation {
namespace method {


using Polygon = util::ConvexSphericalPolygon;
using PolygonArray = std::vector<util::ConvexSphericalPolygon>;


ConservativeSphericalPolygonInterpolationLimiter::
ConservativeSphericalPolygonInterpolationLimiter(const ConservativeSphericalPolygonInterpolation& interpolation):
    interpolation_(interpolation), src_cell_data_(interpolation.src_cell_data_), tgt_cell_data_(interpolation.tgt_cell_data_),
    limiter_(interpolation.limiter_), order_(interpolation.order_),
    src_fs_(interpolation.source()), tgt_fs_(interpolation.target()), data_(interpolation.data_),
    src_points_(data_->src_.points), src_iparam_(data_->src_iparam_), tgt_iparam_(data_->tgt_iparam_),
    tgt_areas_(data_->tgt_.areas) {
    limiter_output_ = interpolation_.limiter_output_;
    limiter_detector_size_ = interpolation_.limiter_detector_size_;
    Log::info() << "limiter output : " << limiter_output_ << std::endl;
    Log::info() << "limiter detector size : " << limiter_detector_size_ << std::endl;
}


double ConservativeSphericalPolygonInterpolationLimiter::limit(const Field& src_field, Field& tgt_field) {
    double mass_change = 0.;
    if (order_ != 2 || limiter_ == "none") {
        return mass_change;
    }
    const auto src_vals = array::make_view<double, 1>(src_field);
    auto tgt_vals       = array::make_view<double, 1>(tgt_field);
    Field tgt_lim_field = tgt_fs_.createField<double>();
    auto tgt_lim_vals   = array::make_view<double, 1>(tgt_lim_field);
    for (idx_t tpt = 0; tpt < tgt_lim_vals.size(); ++tpt) {
        tgt_lim_vals(tpt) = 0.;
    }
    double tgt_smin;
    double tgt_smax;
    if (limiter_ == "clip") {
        for (idx_t tcsp = 0; tcsp < data_->tgt_.csp_size; ++tcsp) {
            idx_t tpt = tgt_cell_data_ ? interpolation_.csp_to_cell(tcsp, data_->tgt_) : data_->tgt_.csp2node[tcsp];
            tgt_smin = std::numeric_limits<double>::max();
            tgt_smax = -tgt_smin;
            const auto& iparam = tgt_iparam_[tcsp];
            if (violation_detected(tpt, iparam, src_vals, tgt_vals, tgt_smin, tgt_smax)) {
                if (limiter_output_ == "points") {
                    tgt_lim_vals(tpt) = 1.;
                    continue;
                }
                if (tgt_vals(tpt) > tgt_smax) {
                    tgt_lim_vals(tpt) = tgt_smax - tgt_vals(tpt);
                }
                else if (tgt_vals(tpt) < tgt_smin) {
                    tgt_lim_vals(tpt) = tgt_smin - tgt_vals(tpt);
                }
            }
        }
    }
    if (limiter_ == "ilmc") {
        const double inf = std::numeric_limits<double>::max();
        std::vector<double> tgt_work_vals(tgt_vals.size());
        std::vector<double> bounds_min(tgt_vals.size(), inf);
        std::vector<double> bounds_max(tgt_vals.size(), -inf);
        std::vector<bool> has_bounds(tgt_vals.size(), false);

        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            tgt_work_vals[tpt] = tgt_vals(tpt);
        }

        for (idx_t tcsp = 0; tcsp < data_->tgt_.csp_size; ++tcsp) {
            idx_t tpt = tgt_cell_data_ ? interpolation_.csp_to_cell(tcsp, data_->tgt_) : data_->tgt_.csp2node[tcsp];
            const auto& iparam = tgt_iparam_[tcsp];
            double smin_loc;
            double smax_loc;
            if (violation_detected(tpt, iparam, src_vals, tgt_vals, smin_loc, smax_loc)) {
                if (has_bounds[tpt]) {
                    bounds_min[tpt] = std::min(bounds_min[tpt], smin_loc);
                    bounds_max[tpt] = std::max(bounds_max[tpt], smax_loc);
                }
                else {
                    has_bounds[tpt] = true;
                    bounds_min[tpt] = smin_loc;
                    bounds_max[tpt] = smax_loc;
                }
            }
        }

        constexpr double eps = std::numeric_limits<double>::epsilon();
        constexpr std::size_t ilmc_max_shells = 3;
        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            if (! has_bounds[tpt]) {
                continue;
            }
            const double bounded_value = std::min(std::max(tgt_work_vals[tpt], bounds_min[tpt]), bounds_max[tpt]);
            double delta_mass = (tgt_work_vals[tpt] - bounded_value) * tgt_areas_[tpt];
            if (std::abs(delta_mass) <= eps) {
                continue;
            }

            tgt_work_vals[tpt] = bounded_value;
            double residual_mass;
            for (int shell = 0; shell < ilmc_max_shells && std::abs(residual_mass) > eps; ++shell) {
                residual_mass = redistribute_local_mass(tpt, delta_mass, bounds_min, bounds_max, has_bounds, tgt_work_vals);
            }
            if (std::abs(residual_mass) > eps) {
                Log::info() << "Residual mass after ILMC redistribution at target point " << tpt << " is " << residual_mass
                          << std::endl;
                tgt_work_vals[tpt] += residual_mass / tgt_areas_[tpt];
            }
        }

        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            if (limiter_output_ == "points") {
                tgt_lim_vals(tpt) = has_bounds[tpt] ? 1. : 0.;
            }
            else {
                tgt_lim_vals(tpt) = tgt_work_vals[tpt] - tgt_vals(tpt);
            }
        }
    }
    if (limiter_ == "zeroslope") {
        scsp_acted_on_tcsp_.clear();
        scsp_acted_on_tcsp_.resize(data_->src_.csp_size);
        compute_src_grad(src_vals);
        std::set<idx_t> send_marked_spt_set;
        for (idx_t tcsp = 0; tcsp < data_->tgt_.csp_size; ++tcsp) {
            idx_t tpt;
            if (tgt_cell_data_) {
                tpt = interpolation_.csp_to_cell(tcsp, data_->tgt_);
            }
            else {
                tpt = data_->tgt_.csp2node[tcsp];
            }
            const auto& iparam = tgt_iparam_[tcsp];
            if (violation_detected(tpt, iparam, src_vals, tgt_vals, tgt_smin, tgt_smax)) {
                for (auto scsp_id : iparam.csp_ids) {
                    idx_t spt;
                    if (src_cell_data_) {
                        spt  = interpolation_.csp_to_cell(scsp_id, data_->src_);
                    }
                    else {
                        spt  = data_->src_.csp2node[scsp_id];
                    }
                    send_marked_spt_set.insert(spt);
                }
                if (limiter_output_ == "points") {
                    tgt_lim_vals(tpt) = 1.;
                }
            }
        }
        auto& mpi_comm = mpi::comm();
        auto mpi_size = mpi_comm.size();
        std::vector<gidx_t> send_marked_spt;
        send_marked_spt.reserve(send_marked_spt_set.size());
        const auto src_global_index = array::make_view<gidx_t, 1>(src_fs_.global_index());
        for (auto idx : send_marked_spt_set) {
            gidx_t gidx = src_global_index(idx);
            send_marked_spt.emplace_back(gidx);
        }
        eckit::mpi::Buffer<gidx_t> recv_marked_scells_buf(mpi_size);
        mpi_comm.allGatherv(send_marked_spt.begin(), send_marked_spt.end(), recv_marked_scells_buf);

        std::unordered_map<gidx_t, idx_t> recv_loc_scells;
        ATLAS_TRACE_SCOPE("Build global-to-local map for source mesh in ConservativeSphericalPolygonInterpolationLimiter") {
            for (idx_t spt_loc = 0; spt_loc < src_global_index.size(); ++spt_loc) {
                auto spt_glo = src_global_index(spt_loc);
                if (recv_loc_scells.find(spt_glo) != recv_loc_scells.end()) {
                    continue;
                }
                recv_loc_scells[spt_glo] = spt_loc;
            }
        }

        auto recv_size = std::accumulate(recv_marked_scells_buf.counts.begin(), recv_marked_scells_buf.counts.end(), 0);
        for (idx_t i_gid = 0; i_gid < recv_size; ++i_gid) {
            idx_t spt = recv_loc_scells[recv_marked_scells_buf.buffer[i_gid]];
            if (src_cell_data_) {
                idx_t scsp_id = spt;  // TODO: convert scell to scsp_id
                limit_contrib_from_source(scsp_id, src_field, tgt_lim_vals);
            }
            else {
                for (auto scsp_id : data_->src_.node2csp[spt]) {
                    limit_contrib_from_source(scsp_id, src_field, tgt_lim_vals);
                }
            }
        }
    }
    if (limiter_output_ == "target") {
        for (idx_t tpt = 0 ; tpt < tgt_vals.size(); ++tpt) {
            mass_change += tgt_lim_vals(tpt)  * tgt_areas_[tpt];
            tgt_vals(tpt) += tgt_lim_vals(tpt);
        }
    }
    else {
        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            tgt_vals(tpt) = tgt_lim_vals(tpt);
        }
    }
    return mass_change;
}


std::vector<idx_t> ConservativeSphericalPolygonInterpolationLimiter::target_neighbours(
    idx_t tpt, ConservativeSphericalPolygonInterpolation::Workspace_get_cell_neighbours& w_cell,
    ConservativeSphericalPolygonInterpolation::Workspace_get_node_neighbours& w_node) const {
    return tgt_cell_data_ ? interpolation_.get_cell_neighbours(interpolation_.tgt_mesh_, tpt, w_cell)
                          : interpolation_.get_node_neighbours(interpolation_.tgt_mesh_, tpt, w_node);
}


double ConservativeSphericalPolygonInterpolationLimiter::redistribute_local_mass(
    idx_t tpt, double delta_mass, const std::vector<double>& smin, const std::vector<double>& smax,
    const std::vector<bool>& has_bounds, std::vector<double>& tgt_work_vals) const {
    const bool distribute_excess = (delta_mass > 0.);
    double remaining = std::abs(delta_mass);
    constexpr double eps = std::numeric_limits<double>::epsilon();
    if (remaining <= eps) {
        return 0.;
    }

    std::vector<bool> visited(tgt_work_vals.size(), false);
    visited[tpt] = true;
    std::vector<idx_t> frontier(1, tpt);
    std::vector<idx_t> next_frontier;

    ConservativeSphericalPolygonInterpolation::Workspace_get_cell_neighbours w_cell;
    ConservativeSphericalPolygonInterpolation::Workspace_get_node_neighbours w_node;

    next_frontier.clear();
    for (auto center : frontier) {
        auto neighbours = target_neighbours(center, w_cell, w_node);
        for (auto nb : neighbours) {
            // if (nb < 0 || nb >= static_cast<idx_t>(tgt_work_vals.size())) {
            //     continue;
            // }
            if (! visited[nb]) {
                visited[nb] = true;
                next_frontier.push_back(nb);
            }
        }
    }
    // if (next_frontier.empty()) {
    //     break;
    // }
    std::vector<std::pair<idx_t, double>> capacities;
    capacities.reserve(next_frontier.size());
    double total_capacity = 0.;
    for (auto nb : next_frontier) {
        if (! has_bounds[nb]) {
            continue;
        }
        double capacity = 0.;
        if (distribute_excess) {
            capacity = std::max((smax[nb] - tgt_work_vals[nb]) * tgt_areas_[nb], 0.);
        }
        else {
            capacity = std::max((tgt_work_vals[nb] - smin[nb]) * tgt_areas_[nb], 0.);
        }
        if (capacity > 0.) {
            capacities.emplace_back(nb, capacity);
            total_capacity += capacity;
        }
    }
    if (total_capacity > 0.) {
        const double transferred_mass = std::min(remaining, total_capacity);
        for (const auto& cap : capacities) {
            idx_t nb = cap.first;
            const double dm = transferred_mass * cap.second / total_capacity;
            if (distribute_excess) {
                tgt_work_vals[nb] += dm / tgt_areas_[nb];
            }
            else {
                tgt_work_vals[nb] -= dm / tgt_areas_[nb];
            }
        }
        remaining -= transferred_mass;
    }
    frontier.swap(next_frontier);

    if (remaining > eps) {
        std::vector<std::pair<idx_t, double>> capacities;
        capacities.reserve(tgt_work_vals.size());
        double total_capacity = 0.;
        for (idx_t idx = 0; idx < static_cast<idx_t>(tgt_work_vals.size()); ++idx) {
            if (idx == tpt || !has_bounds[idx] || tgt_areas_[idx] <= 0.) {
                continue;
            }
            double capacity = 0.;
            if (distribute_excess) {
                capacity = std::max((smax[idx] - tgt_work_vals[idx]) * tgt_areas_[idx], 0.);
            }
            else {
                capacity = std::max((tgt_work_vals[idx] - smin[idx]) * tgt_areas_[idx], 0.);
            }
            if (capacity > 0.) {
                capacities.emplace_back(idx, capacity);
                total_capacity += capacity;
            }
        }

        if (total_capacity > 0.) {
            const double transferred_mass = std::min(remaining, total_capacity);
            for (const auto& cap : capacities) {
                idx_t idx = cap.first;
                const double dm = transferred_mass * cap.second / total_capacity;
                if (distribute_excess) {
                    tgt_work_vals[idx] += dm / tgt_areas_[idx];
                }
                else {
                    tgt_work_vals[idx] -= dm / tgt_areas_[idx];
                }
            }
            remaining -= transferred_mass;
        }
    }

    return distribute_excess ? remaining : -remaining;
}


bool ConservativeSphericalPolygonInterpolationLimiter::
violation_detected(idx_t tpt, const InterpolationParameters& tiparam, const array::ArrayView<double,1>& src_vals,
    const array::ArrayView<double,1>& tgt_vals, double& smin, double& smax) const {
    smin = std::numeric_limits<double>::max();
    smax = -smin;
    ConservativeSphericalPolygonInterpolation::Workspace_get_cell_neighbours w_cell; // TODO: move even higher up in scope ??
    ConservativeSphericalPolygonInterpolation::Workspace_get_node_neighbours w_node;
    bool do_neighbours = (tiparam.csp_ids.size() <= limiter_detector_size_);
    for (idx_t i_scsp = 0; i_scsp < tiparam.csp_ids.size(); ++i_scsp) {
        idx_t scsp_id = tiparam.csp_ids[i_scsp];
        idx_t spt;
        std::vector<idx_t> src_neighbours;
        spt  = src_cell_data_ ? interpolation_.csp_to_cell(scsp_id, data_->src_) : data_->src_.csp2node[scsp_id];
        smax = std::max(smax, src_vals(spt));
        smin = std::min(smin, src_vals(spt));
        if (do_neighbours) {
            src_neighbours = src_cell_data_ ? interpolation_.get_cell_neighbours(interpolation_.src_mesh_, spt, w_cell) :
                interpolation_.get_node_neighbours(interpolation_.src_mesh_, spt, w_node);
            for (auto nb : src_neighbours) {
                smax = std::max(smax, src_vals(nb));
                smin = std::min(smin, src_vals(nb));
            }
        }
    }
    if (smin == std::numeric_limits<double>::max()) {
        // do not trigger limiting for ghost and halo target points which are not associated with any source polygon
        return false;
    }
    bool undershoot = (tgt_vals(tpt) < smin);
    bool overshoot = (tgt_vals(tpt) > smax);
    return (undershoot || overshoot);
}


void ConservativeSphericalPolygonInterpolationLimiter::
compute_src_grad(const array::ArrayView<double,1>& src_vals) {
    ATLAS_TRACE("Compute source gradients");
    src_grads_.resize(src_vals.size());
    if (src_cell_data_) {
        for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
            src_grads_[scell] = interpolation_.src_gradient_celldata(scell, src_vals);
        }
    }
    else {
        for (idx_t snode = 0; snode < src_vals.size(); ++snode) {
                src_grads_[snode] = interpolation_.src_gradient_nodedata(snode, src_vals);
        }
    }
}


void ConservativeSphericalPolygonInterpolationLimiter::
limit_contrib_from_source(idx_t scsp_id, const Field& src_field, array::ArrayView<double,1>& tgt_lim_vals) {
    idx_t spt;
    if (src_cell_data_) {
        spt = interpolation_.csp_to_cell(scsp_id, data_->src_);
    }
    else {
        spt = data_->src_.csp2node[scsp_id];
    }
    const PointXYZ& src_barycentre = src_points_[spt];
    PointXYZ spt_grad  = src_grads_[spt];
    spt_grad           = spt_grad - PointXYZ::mul(src_barycentre, PointXYZ::dot(spt_grad, src_barycentre));
    auto& siparam = src_iparam_[scsp_id];
    for (idx_t i_tcsp_collateral = 0; i_tcsp_collateral < siparam.csp_ids.size(); ++i_tcsp_collateral) {
        auto tcsp_collateral = siparam.csp_ids[i_tcsp_collateral];
        const auto& iparam_collateral = tgt_iparam_[tcsp_collateral];
        idx_t tpt_collateral;
        if (tgt_cell_data_) {
            tpt_collateral = interpolation_.csp_to_cell(tcsp_collateral, data_->tgt_);
        }
        else {
            tpt_collateral = data_->tgt_.csp2node[tcsp_collateral];
        }
        // find the index of the scsp_id entry in iparam_collateral.csp_ids
        auto scsp_it = std::find(iparam_collateral.csp_ids.begin(), iparam_collateral.csp_ids.end(), scsp_id);
        ATLAS_ASSERT(scsp_it != iparam_collateral.csp_ids.end());
        idx_t scsp_idx = scsp_it - iparam_collateral.csp_ids.begin();
        ATLAS_ASSERT(iparam_collateral.csp_ids[scsp_idx] == scsp_id);
        double tgt_lim_val = iparam_collateral.weights[scsp_idx] * PointXYZ::dot(spt_grad, iparam_collateral.centroids[scsp_idx] - src_barycentre);
        if (tgt_areas_[tpt_collateral] > 0.) {
            tgt_lim_val /= tgt_areas_[tpt_collateral];
        }
        auto& it = scsp_acted_on_tcsp_[scsp_id];
        if (std::find(it.tcsp_done.begin(), it.tcsp_done.end(), tcsp_collateral) == it.tcsp_done.end()) {
            it.tcsp_done.push_back(tcsp_collateral);
            if (limiter_output_ == "points") {
                if (tgt_lim_vals(tpt_collateral) < 0.5) {
                    tgt_lim_vals(tpt_collateral) = -1.;
                }
            }
            else {
                tgt_lim_vals(tpt_collateral) -= tgt_lim_val;
            }
        }
    }
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas