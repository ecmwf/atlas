/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include "atlas/functionspace.h"
#include "atlas/interpolation/method/Method.h"
#include "atlas/interpolation/method/unstructured/ConservativeSphericalPolygonInterpolation.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/ConvexSphericalPolygon.h"


namespace atlas {
namespace interpolation {
namespace method {


using Indices = std::vector<idx_t>;


class ConservativeSphericalPolygonInterpolationLimiter {
public:
    using Data = ConservativeSphericalPolygonInterpolation::Data;
    using InterpolationParameters = ConservativeSphericalPolygonInterpolation::InterpolationParameters;
    using Polygon = util::ConvexSphericalPolygon;
    using PolygonArray = std::vector<util::ConvexSphericalPolygon>;
    struct SrcActed {
        Indices tcsp_done;
    };

    ConservativeSphericalPolygonInterpolationLimiter(const ConservativeSphericalPolygonInterpolation& interpolation);

    const ConservativeSphericalPolygonInterpolation& interpolation() const { return interpolation_; }
    double limit(const Field& src_field, Field& tgt_field);

    // interpolation::Cache createCache() const override;

private:
    bool violation_detected(idx_t tcell, const InterpolationParameters& tiparam, const array::ArrayView<double,1>& src_vals,
        const array::ArrayView<double,1>& tgt_vals, double& smin, double& smax) const;
    void compute_src_grad(const array::ArrayView<double,1>& src_vals);
    void limit_contrib_from_source(idx_t scsp_id, const Field& src_field, array::ArrayView<double,1>& tgt_lim_vals);

private:
    const ConservativeSphericalPolygonInterpolation& interpolation_;
    bool src_cell_data_;
    bool tgt_cell_data_;
    std::vector<PointXYZ> src_grads_;
    std::vector<SrcActed> spt_acted_on_tcsp_;
    std::string limiter_;
    int order_;
    const FunctionSpace src_fs_;
    const FunctionSpace tgt_fs_;
    unsigned int limiter_override_tgt_ = 0;

    // Cache cache_;                          // Storage of cache if any was passed to constructor
    // std::shared_ptr<Data> sharable_data_;  // Storage of new data_, only allocated if cache is empty
    const ConservativeSphericalPolygonInterpolation::Data* data_;  // Read-only access to data, pointing either to cache_ or sharable_data_
    const std::vector<PointXYZ>& src_points_;
    const std::vector<InterpolationParameters>& src_iparam_;
    const std::vector<InterpolationParameters>& tgt_iparam_;
    const std::vector<double>& tgt_areas_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
