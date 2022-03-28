/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/knn/GridBoxMaximum.h"

#include <algorithm>
#include <limits>
#include <vector>

#include "eckit/log/ProgressTimer.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/runtime/Exception.h"


namespace atlas {
namespace interpolation {
namespace method {


namespace {
MethodBuilder<GridBoxMaximum> __builder("grid-box-maximum");
}


void GridBoxMaximum::do_execute(const FieldSet& source, FieldSet& target, Metadata& metadata) const {
    ATLAS_ASSERT(source.size() == target.size());

    for (idx_t i = 0; i < source.size(); ++i) {
        do_execute(source[i], target[i], metadata);
    }
}


void GridBoxMaximum::do_execute(const Field& source, Field& target, Metadata&) const {
    ATLAS_TRACE("atlas::interpolation::method::GridBoxMaximum::do_execute()");

    // set arrays
    ATLAS_ASSERT(source.rank() == 1);
    ATLAS_ASSERT(target.rank() == 1);

    auto xarray = atlas::array::make_view<double, 1>(source);
    auto yarray = atlas::array::make_view<double, 1>(target);
    ATLAS_ASSERT(xarray.size() == sourceBoxes_.size());
    ATLAS_ASSERT(yarray.size() == targetBoxes_.size());

    yarray.assign(0.);
    failures_.clear();


    if (!matrixFree_) {
        const Matrix& m = matrix();
        Matrix::const_iterator k(m);

        for (decltype(m.rows()) i = 0, j = 0; i < m.rows(); ++i) {
            double max = std::numeric_limits<double>::lowest();
            bool found = false;

            for (; k != m.end(i); ++k) {
                ATLAS_ASSERT(k.col() < size_t(xarray.shape(0)));
                auto value = xarray[k.col()];
                if (max < value) {
                    max = value;
                    j   = k.col();
                }
                found = true;
            }

            ATLAS_ASSERT(found);
            yarray[i] = xarray[j];
        }
        return;
    }


    // ensure GridBoxMethod::setup()
    functionspace::PointCloud tgt = target_;
    ATLAS_ASSERT(tgt);

    ATLAS_ASSERT(searchRadius_ > 0.);
    ATLAS_ASSERT(!sourceBoxes_.empty());
    ATLAS_ASSERT(!targetBoxes_.empty());


    // interpolate
    eckit::ProgressTimer progress("Intersecting", targetBoxes_.size(), "grid box", double(5.));

    std::vector<Triplet> triplets;
    size_t i = 0;
    for (auto p : tgt.iterate().lonlat()) {
        ++progress;

        if (intersect(i, targetBoxes_.at(i), pTree_.closestPointsWithinRadius(p, searchRadius_), triplets)) {
            auto triplet = std::max_element(triplets.begin(), triplets.end(), [](const Triplet& a, const Triplet& b) {
                return !eckit::types::is_approximately_greater_or_equal(a.value(), b.value());
            });

            yarray[i] = xarray[triplet->col()];
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
