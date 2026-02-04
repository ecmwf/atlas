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
#include "atlas/util/ConvexSphericalPolygon.h"


namespace atlas {
namespace interpolation {
namespace method {


using Indices = std::vector<idx_t>;


class ConservativeSphericalPolygonInterpolationLimiter {
public:
    using Data = ConservativeSphericalPolygonInterpolation::Data;
    using InterpolationParameters = ConservativeSphericalPolygonInterpolation::InterpolationParameters;

public:
    ConservativeSphericalPolygonInterpolationLimiter(const ConservativeSphericalPolygonInterpolation& interpolation);

    const ConservativeSphericalPolygonInterpolation& interpolation() const { return interpolation_; }
    void limit(const Field& src_field, Field& tgt_field);

    // interpolation::Cache createCache() const override;

private:

    using Polygon = util::ConvexSphericalPolygon;
    using PolygonArray = std::vector<util::ConvexSphericalPolygon>;

private:
    const ConservativeSphericalPolygonInterpolation& interpolation_;
    bool src_cell_data_;
    bool tgt_cell_data_;
    std::string limiter_;
    int order_;
    bool matrix_free_;
    const FunctionSpace src_fs_;
    const FunctionSpace tgt_fs_;

    // Cache cache_;                          // Storage of cache if any was passed to constructor
    // std::shared_ptr<Data> sharable_data_;  // Storage of new data_, only allocated if cache is empty
    const ConservativeSphericalPolygonInterpolation::Data* data_;  // Read-only access to data, pointing either to cache_ or sharable_data_
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
