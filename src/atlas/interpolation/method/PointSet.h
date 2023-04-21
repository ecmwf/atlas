/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <stdint.h>
#include <cassert>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "eckit/config/Resource.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/interpolation/method/PointIndex3.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace interpolation {
namespace method {

//----------------------------------------------------------------------------------------------------------------------

class PointSet {
public:  // types
    typedef PointIndex3::Point Point;
    typedef PointIndex3::iterator iterator;
    typedef std::map<size_t, size_t> DupStore_t;

public:  // methods
    PointSet(const std::vector<Point>& ipts);

    PointSet(atlas::Mesh& mesh);

    ~PointSet() { delete tree_; }

    DupStore_t& duplicates() { return duplicates_; }

    template <typename POINT_T>
    void list_unique_points(std::vector<POINT_T>& opts) {
        ATLAS_TRACE("Finding unique points");

        ATLAS_ASSERT(opts.empty());

        opts.reserve(npts_);

        for (PointIndex3::iterator i = tree_->begin(); i != tree_->end(); ++i) {
            Point p(i->point());
            size_t ip = i->payload();
            //            std::cout << "point " << ip << " " << p << std::endl;
            size_t uidx = unique(p, ip);
            if (ip == uidx) {
                opts.push_back(POINT_T(p.data()));
                //                std::cout << "----> UNIQ " << ip << std::endl;
            }
            else {
                //                std::cout << "----> DUP " << ip << " -> " << uidx <<
                //                std::endl;
            }
            //            ++show_progress;
        }
    }

    void list_unique_points(std::vector<size_t>& opts) {
        ATLAS_TRACE("Finding unique points");

        ATLAS_ASSERT(opts.empty());

        opts.reserve(npts_);

        for (PointIndex3::iterator i = tree_->begin(); i != tree_->end(); ++i) {
            Point p(i->point());
            size_t ip = i->payload();
            //            std::cout << "point " << ip << " " << p << std::endl;
            size_t uidx = unique(p, ip);
            if (ip == uidx) {
                opts.push_back(ip);
                //                std::cout << "----> UNIQ " << ip << std::endl;
            }
            else {
                //                std::cout << "----> DUP " << ip << " -> " << uidx <<
                //                std::endl;
            }
            //            ++show_progress;
        }
    }


    size_t unique(const Point& p, size_t idx = std::numeric_limits<size_t>::max()) {
        DupStore_t::iterator dit = duplicates_.find(idx);
        if (dit != duplicates_.end()) {
            //                std::cout << "      !! DUPLICATE !!" << std::endl;
            return dit->second;
        }
        else
            return this->search_unique(p, idx, 0);
    }

    iterator begin() { return tree_->begin(); }
    iterator end() { return tree_->end(); }

    size_t size() const { return npts_; }

protected:  // methods
    template <typename V>
    void build(const V& ipts) {
        static bool fastBuildKDTrees = eckit::Resource<bool>("$ATLAS_FAST_BUILD_KDTREES", true);
        tree_                        = new PointIndex3();

        if (fastBuildKDTrees) {
            std::vector<PointIndex3::Value> pidx;
            pidx.reserve(npts_);

            for (size_t ip = 0; ip < npts_; ++ip) {
                pidx.push_back(PointIndex3::Value(PointIndex3::Point(ipts[ip]), ip));
            }

            tree_->build(pidx.begin(), pidx.end());
        }
        else {
            for (size_t ip = 0; ip < npts_; ++ip) {
                tree_->insert(PointIndex3::Value(PointIndex3::Point(ipts[ip]), ip));
            }
        }
    }

    size_t search_unique(const Point& p, size_t idx, uint32_t n);

protected:
    size_t Kn(uint32_t n) {
        if (!n)
            return 2;
        else
            return 2 + std::pow(2, n) * 180;
    }

    bool duplicate(size_t idx) { return duplicates_.find(idx) != duplicates_.end(); }

private:
    size_t npts_;
    PointIndex3* tree_;
    DupStore_t duplicates_;  ///< map from duplicate idx to idx representing group
                             /// of points
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
