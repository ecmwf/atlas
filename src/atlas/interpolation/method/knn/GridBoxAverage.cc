/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/knn/GridBoxAverage.h"

#include <vector>

#include "eckit/log/ProgressTimer.h"

#include "atlas/array.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/runtime/Exception.h"


namespace atlas {
namespace interpolation {
namespace method {


namespace {
MethodBuilder<GridBoxAverage> __builder("grid-box-average");
}


void GridBoxAverage::do_execute(const FieldSet& source, FieldSet& target, Metadata& metadata) const {
    ATLAS_ASSERT(source.size() == target.size());

    // Matrix-based interpolation is handled by base (Method) class
    // TODO: exploit sparse/dense matrix multiplication
    for (idx_t i = 0; i < source.size(); ++i) {
        if (matrixFree_) {
            GridBoxAverage::do_execute(source[i], target[i], metadata);
        }
        else {
            Method::do_execute(source[i], target[i], metadata);
        }
    }
}


void GridBoxAverage::do_execute(const Field& source, Field& target, Metadata& metadata) const {
    ATLAS_TRACE("atlas::interpolation::method::GridBoxAverage::do_execute()");

    // Matrix-based interpolation is handled by base (Method) class
    if (!matrixFree_) {
        Method::do_execute(source, target, metadata);
        return;
    }


    // ensure GridBoxMethod::setup()
    functionspace::PointCloud tgt = target_;
    ATLAS_ASSERT(tgt);

    ATLAS_ASSERT(searchRadius_ > 0.);
    ATLAS_ASSERT(!sourceBoxes_.empty());
    ATLAS_ASSERT(!targetBoxes_.empty());


    // set arrays
    ATLAS_ASSERT(source.rank() == 1);
    ATLAS_ASSERT(target.rank() == 1);

    auto xarray = atlas::array::make_view<double, 1>(source);
    auto yarray = atlas::array::make_view<double, 1>(target);
    ATLAS_ASSERT(xarray.size() == sourceBoxes_.size());
    ATLAS_ASSERT(yarray.size() == targetBoxes_.size());

    yarray.assign(0.);
    failures_.clear();


    // interpolate
    eckit::ProgressTimer progress("Intersecting", targetBoxes_.size(), "grid box", double(5.));

    std::vector<Triplet> triplets;
    size_t i = 0;
    for (auto p : tgt.iterate().lonlat()) {
        ++progress;

        if (intersect(i, targetBoxes_.at(i), pTree_.closestPointsWithinRadius(p, searchRadius_), triplets)) {
            auto& y = yarray[i];
            for (auto& t : triplets) {
                y += xarray[t.col()] * t.value();
            }
        }

        ++i;
    }

    if (!failures_.empty()) {
        giveUp(failures_);
    }
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
