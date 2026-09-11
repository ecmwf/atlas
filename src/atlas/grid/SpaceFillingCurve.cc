/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/SpaceFillingCurve.h"

#include <cmath>
#include <limits>

namespace atlas {
namespace grid {

// -------------------------------------------------------------------------------------

HilbertCurve::HilbertCurve(const Domain& domain, idx_t levels): domain_{domain}, max_level_(levels) {
    nb_keys_2_ = gidx_t(std::pow(gidx_t(4), gidx_t(max_level_)));
    nb_keys_   = nb_keys_2_ * 2;
}

gidx_t HilbertCurve::operator()(const PointXY& point) const {
    box_t box;
    box[A]            = {domain_.xmin(), domain_.ymax()};
    box[B]            = {domain_.xmin(), domain_.ymin()};
    box[C]            = {domain_.xmax(), domain_.ymin()};
    box[D]            = {domain_.xmax(), domain_.ymax()};
    const double xmid = (domain_.xmin() + domain_.xmax()) * 0.5;
    if (point.x() < xmid) {
        box[C].x() = xmid;
        box[D].x() = xmid;
        return recursive_algorithm(point, box, 0);
    }
    else {
        box[A].x() = xmid;
        box[B].x() = xmid;
        return recursive_algorithm(point, box, 0) + nb_keys_2_;
    }
}

gidx_t HilbertCurve::recursive_algorithm(const PointXY& p, const box_t& box, idx_t level) const {
    if (level == max_level_) {
        return 0;
    }

    double min_distance = std::numeric_limits<double>::max();

    auto compute_distance2 = [](const PointXY& p1, const PointXY& p2) {
        // workaround because of eckit 1.3.2 issue with constness in KPoint
        double d = 0;
        for (size_t i = 0; i < 2; i++) {
            double dx = p1[i] - p2[i];
            d += dx * dx;
        }
        return d;
    };

    auto compute_average = [](const PointXY& p1, const PointXY& p2) {
        // workaround because of eckit 1.3.2 issue with constness in KPoint
        PointXY avg;
        avg.x() = p1.x() + p2.x();
        avg.x() *= 0.5;
        avg.y() = p1.y() + p2.y();
        avg.y() *= 0.5;
        return avg;
    };

    idx_t quadrant{0};
    for (idx_t idx = 0; idx < 4; ++idx) {
        // double distance = box[idx].distance2( p );  // does not compile with eckit 1.3.2
        double distance = compute_distance2(p, box[idx]);  // workaround
        if (distance < min_distance) {
            quadrant     = idx;
            min_distance = distance;
        }
    }

    box_t box_quadrant;
    switch (quadrant) {
        case A:
            box_quadrant[A] = box[A];
            box_quadrant[B] = compute_average(box[A], box[D]);  // workaround
            box_quadrant[C] = compute_average(box[A], box[C]);  // workaround
            box_quadrant[D] = compute_average(box[A], box[B]);  // workaround
            break;
        case B:
            box_quadrant[B] = box[B];
            box_quadrant[A] = compute_average(box[B], box[A]);  // workaround
            box_quadrant[C] = compute_average(box[B], box[C]);  // workaround
            box_quadrant[D] = compute_average(box[B], box[D]);  // workaround
            break;
        case C:
            box_quadrant[C] = box[C];
            box_quadrant[A] = compute_average(box[C], box[A]);  // workaround
            box_quadrant[B] = compute_average(box[C], box[B]);  // workaround
            box_quadrant[D] = compute_average(box[C], box[D]);  // workaround
            break;
        case D:
            box_quadrant[D] = box[D];
            box_quadrant[A] = compute_average(box[D], box[C]);  // workaround
            box_quadrant[B] = compute_average(box[D], box[B]);  // workaround
            box_quadrant[C] = compute_average(box[D], box[A]);  // workaround
            break;
    }

    // The key has 4 possible values per recursion (1 for each quadrant),
    // which can be represented by 2 bits per recursion
    //   A --> 00
    //   B --> 01
    //   C --> 10
    //   D --> 11
    // Trailing zero-bits are added depending on the level:
    //   level max_level_-1 --> none
    //   level max_level_-2 --> 00
    //   level max_level_-2 --> 0000
    //   level max_level_-3 --> 000000
    gidx_t key = 0;
    auto index = (max_level_ - level) * 2 - 1;
    gidx_t mask;

    // Create a mask value with all trailing bits for leftmost bit (of 2)
    mask = gidx_t(1) << index;

    // Add mask to key
    if (quadrant == C || quadrant == D) {
        key |= mask;
    }

    // Create a mask value with all trailing bits for rightmost bit (of 2)
    mask = gidx_t(1) << (index - 1);

    // Add mask to key
    if (quadrant == B || quadrant == D) {
        key |= mask;
    }

    return recursive_algorithm(p, box_quadrant, level + 1) + key;
}

// -------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
