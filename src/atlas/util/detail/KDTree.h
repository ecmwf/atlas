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

#include <iosfwd>
#include <memory>

#include "eckit/container/KDTree.h"

#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/Object.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace util {
namespace detail {

//------------------------------------------------------------------------------------------------------

// Abstract KDTree, intended to erase the internal KDTree type from eckit (Mapped or Memory)
// For usage, see atlas::util::KDTree
template <typename PayloadT, typename PointT = Point3>
class KDTreeBase : public Object {
#define ENABLE_IF_3D_AND_IS_LONLAT(POINT) \
    typename = typename std::enable_if < PointT::DIMS == 3 && POINT::DIMS == 2 > ::type

private:
    Geometry geometry_;

public:
    using Payload     = PayloadT;
    using Point       = PointT;
    using PayloadList = std::vector<Payload>;

    struct KDTreeTraits {
        using Point   = PointT;
        using Payload = PayloadT;
    };

    struct Value {
        using Point   = typename KDTreeTraits::Point;
        using Payload = typename KDTreeTraits::Payload;

        template <typename Node>
        Value(const Node& node): Value(node.point(), node.payload(), node.distance()) {}

        Value(const Point& point, const Payload& payload): point_(point), payload_(payload) {}

        Value(const Point& point, const Payload& payload, double distance):
            point_(point), payload_(payload), distance_(distance) {}

        const Point& point() const { return point_; }
        const Payload& payload() const { return payload_; }
        const double& distance() const { return distance_; }

    private:
        Point point_;
        Payload payload_;
        double distance_;
    };

    class ValueList : public std::vector<Value> {
    public:
        PayloadList payloads() const {
            PayloadList list;
            list.reserve(this->size());
            for (auto& item : *this) {
                list.emplace_back(item.payload());
            }
            return list;
        }

        template <typename NodeList>
        ValueList(const NodeList& _list) {
            this->reserve(_list.size());
            for (const auto& item : _list) {
                this->emplace_back(item);
            }
        }

        void print(std::ostream& out) const {
            out << "[";
            for (size_t i = 0; i < this->size(); ++i) {
                const auto& value = this->at(i);
                out << "{payload:" << value.payload() << ",distance:" << value.distance() << ",point:" << value.point()
                    << "}";
                if (i < this->size() - 1) {
                    out << ",";
                }
            }
            out << "]";
        }
        friend std::ostream& operator<<(std::ostream& out, ValueList& This) {
            This.print(out);
            return out;
        }
    };


public:
    KDTreeBase() = default;

    KDTreeBase(const Geometry& geometry): geometry_(geometry) {}

    virtual ~KDTreeBase() = default;

    const Geometry& geometry() const { return geometry_; }

    bool empty() const { return size() == 0; }

    virtual idx_t size() const = 0;

    virtual size_t footprint() const = 0;

    /// @brief Reserve memory for building the kdtree in one shot (optional, at cost of extra memory)
    /// Implementation depends in derived classes
    virtual void reserve(idx_t /*size*/) {}

    /// @brief Insert spherical point (lon,lat) or 3D cartesian point (x,y,z)
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    template <typename Point>
    void insert(const Point& p, const Payload& payload) {
        do_insert(p, payload);
    }

    /// @brief Insert Value
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    virtual void insert(const Value&) = 0;

    /// @brief Build the kd-tree in one shot, if memory has been reserved, depending on derived class implementation
    /// This will need to be called before all search functions like closestPoints().
    virtual void build() {}

    /// @brief Build the kd-tree in one shot
    virtual void build(std::vector<Value>& values) = 0;

    /// @brief Build with spherical points (lon,lat) where longitudes, latitudes, and payloads are separate containers.
    /// Memory will be reserved with reserve() to match the size
    template <typename Longitudes, typename Latitudes, typename Payloads>
    void build(const Longitudes& longitudes, const Latitudes& latitudes, const Payloads& payloads) {
        build(longitudes.begin(), longitudes.end(), latitudes.begin(), latitudes.end(), payloads.begin(),
              payloads.end());
    }

    /// @brief Build with spherical points (lon,lat) given separate iterator ranges for longitudes, latitudes, and payloads.
    /// Memory will be reserved with reserve() to match the size
    template <typename LongitudesIterator, typename LatitudesIterator, typename PayloadsIterator>
    void build(const LongitudesIterator& longitudes_begin, const LongitudesIterator& longitudes_end,
               const LatitudesIterator& latitudes_begin, const LatitudesIterator& latitudes_end,
               const PayloadsIterator& payloads_begin, const PayloadsIterator& payloads_end) {
        reserve(std::distance(payloads_begin, payloads_end));
        auto lon     = longitudes_begin;
        auto lat     = latitudes_begin;
        auto payload = payloads_begin;
        for (; lon != longitudes_end && lat != latitudes_end && payload != payloads_end; lon++, lat++, payload++) {
            insert(Point2{*lon, *lat}, *payload);
        }
        ATLAS_ASSERT(lon == longitudes_end);
        ATLAS_ASSERT(lat == latitudes_end);
        ATLAS_ASSERT(payload == payloads_end);
        build();
    }

    /// @brief Build with spherical points (lon,lat) where longitudes, latitudes, and payloads are separate containers.
    /// Memory will be reserved with reserve() to match the size
    template <typename Points, typename Payloads>
    void build(const Points& points, const Payloads& payloads) {
        build(points.begin(), points.end(), payloads.begin(), payloads.end());
    }

    /// @brief Build with spherical points (lon,lat) given separate iterator ranges for longitudes, latitudes, and payloads.
    /// Memory will be reserved with reserve() to match the size
    template <typename PointIterator, typename PayloadsIterator>
    void build(const PointIterator& points_begin, const PointIterator& points_end,
               const PayloadsIterator& payloads_begin, const PayloadsIterator& payloads_end) {
        reserve(std::distance(points_begin, points_end));
        PointIterator point      = points_begin;
        PayloadsIterator payload = payloads_begin;
        for (; point != points_end && payload != payloads_end; ++point, ++payload) {
            insert(*point, *payload);
        }
        ATLAS_ASSERT(point == points_end);
        ATLAS_ASSERT(payload == payloads_end);
        build();
    }


    /// @brief Find k nearest neighbours given a 3D cartesian point (x,y,z) or 2D lonlat point(lon,lat)
    template <typename Point>
    ValueList closestPoints(const Point& p, size_t k) const {
        return do_closestPoints(p, k);
    }

    /// @brief Find nearest neighbour given a 3D cartesian point (x,y,z)
    template <typename Point>
    Value closestPoint(const Point& p) const {
        return do_closestPoint(p);
    }

    /// @brief Find all points within a distance of given radius from a given point (x,y,z)
    template <typename Point>
    ValueList closestPointsWithinRadius(const Point& p, double radius) const {
        return do_closestPointsWithinRadius(p, radius);
    }

private:
    /// @brief Insert spherical point (lon,lat)
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    void do_insert(const Point& p, const Payload& payload) { insert(Value{p, payload}); }


    /// @brief Insert spherical point (lon,lat)
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    template <typename LonLat, ENABLE_IF_3D_AND_IS_LONLAT(LonLat)>
    void do_insert(const LonLat& p, const Payload& payload) {
        do_insert(make_Point(p), payload);
    }

    /// @brief Find k nearest neighbours given a 3D cartesian point (x,y,z)
    virtual ValueList do_closestPoints(const Point&, size_t k) const = 0;

    /// @brief Find nearest neighbour given a 3D cartesian point (x,y,z)
    virtual Value do_closestPoint(const Point&) const = 0;

    /// @brief Find all points within a distance of given radius from a given point (x,y,z)
    virtual ValueList do_closestPointsWithinRadius(const Point&, double radius) const = 0;


    /// @brief Find k nearest neighbour given a 2D lonlat point (lon,lat)
    template <typename LonLat, ENABLE_IF_3D_AND_IS_LONLAT(LonLat)>
    ValueList do_closestPoints(const LonLat& p, size_t k) const {
        return do_closestPoints(make_Point(p), k);
    }

    /// @brief Find nearest neighbour given a 2D lonlat point (lon,lat)
    template <typename LonLat, ENABLE_IF_3D_AND_IS_LONLAT(LonLat)>
    Value do_closestPoint(const LonLat& p) const {
        return do_closestPoint(make_Point(p));
    }

    /// @brief Find all points within a distance of given radius from a given point (lon,lat)
    template <typename LonLat, ENABLE_IF_3D_AND_IS_LONLAT(LonLat)>
    ValueList do_closestPointsWithinRadius(const LonLat& p, double radius) const {
        return do_closestPointsWithinRadius(make_Point(p), radius);
    }

    template <typename LonLat, ENABLE_IF_3D_AND_IS_LONLAT(LonLat)>
    Point make_Point(const LonLat& lonlat) const {
        static_assert(std::is_base_of<Point2, LonLat>::value, "LonLat must be derived from Point2");
        static_assert(std::is_base_of<Point3, Point>::value, "Point must be derived from Point3");
        Point xyz;
        geometry().lonlat2xyz(lonlat, xyz);
        return xyz;
    }
#undef ENABLE_IF_3D_AND_IS_LONLAT
};

//------------------------------------------------------------------------------------------------------
// Concrete implementation

template <typename TreeT, typename PayloadT = typename TreeT::Payload, typename PointT = typename TreeT::Point>
class KDTree_eckit : public KDTreeBase<PayloadT, PointT> {
    using Tree = TreeT;
    using Base = KDTreeBase<PayloadT, PointT>;

public:
    using Point       = typename Base::Point;
    using Payload     = typename Base::Payload;
    using PayloadList = typename Base::PayloadList;
    using Value       = typename Base::Value;
    using ValueList   = typename Base::ValueList;

    using Base::build;
    using Base::closestPoint;
    using Base::closestPoints;
    using Base::closestPointsWithinRadius;
    using Base::insert;
    using Base::reserve;

public:
    template <typename... Args>
    KDTree_eckit(const Geometry& geometry, Args... args):
        KDTreeBase<Payload, Point>(geometry), tree_(new Tree(args...)) {
        static_asserts();
    }


    template <typename... Args>
    KDTree_eckit(Args... args): tree_(new Tree(args...)) {
        static_asserts();
    }

    KDTree_eckit(const std::shared_ptr<Tree>& tree): tree_(tree) { static_asserts(); }

    KDTree_eckit(const std::shared_ptr<Tree>& tree, const Geometry& geometry):
        KDTreeBase<Payload, Point>(geometry), tree_(tree) {
        static_asserts();
    }

    idx_t size() const override {
#if ATLAS_ECKIT_VERSION_AT_LEAST(1, 13, 2)
        return static_cast<idx_t>(tree_->size());
#else
        // Assume ECKIT-515 not implemented.
        // Very bad implementation, with a workaround for empty tree
        idx_t size{0};
        try {
            for (const auto& item : *tree_) {
                size++;
            }
        }
        catch (const eckit::AssertionFailed&) {
        }
        return size;
#endif
    }

    size_t footprint() const override { return size() * sizeof(typename Tree::Node); }

    void reserve(idx_t size) override;

    void build() override;

    void build(std::vector<Value>&) override;

    /// @brief Insert 3D cartesian point (x,y,z)
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    void insert(const Value& value) override;

    /// @brief Find k nearest neighbours given a 3D cartesian point (x,y,z)
    ValueList do_closestPoints(const Point&, size_t k) const override;

    /// @brief Find nearest neighbour given a 3D cartesian point (x,y,z)
    Value do_closestPoint(const Point&) const override;

    /// @brief Find all points within a distance of given radius from a given point (x,y,z)
    ValueList do_closestPointsWithinRadius(const Point&, double radius) const override;

    const Tree& tree() const { return *tree_; }

private:
    void assert_built() const;

    static void static_asserts() {
        static_assert(std::is_convertible<typename Tree::Payload, Payload>::value,
                      "Tree::Payload must be convertible this Payload type");
        static_assert(std::is_convertible<typename Tree::Point, Point>::value,
                      "Tree::Point must be convertible to this Point type");
    }

private:
    std::vector<Value> tmp_;
    mutable std::shared_ptr<Tree> tree_;  // mutable because its member functions are non-const...
};

//------------------------------------------------------------------------------------------------------

template <typename Payload, typename Point>
using KDTreeMemory = KDTree_eckit<typename eckit::KDTreeMemory<typename KDTreeBase<Payload, Point>::KDTreeTraits>>;

//------------------------------------------------------------------------------------------------------

template <typename TreeT, typename PayloadT, typename PointT>
void KDTree_eckit<TreeT, PayloadT, PointT>::reserve(idx_t size) {
    tmp_.reserve(size);
}


template <typename TreeT, typename PayloadT, typename PointT>
void KDTree_eckit<TreeT, PayloadT, PointT>::insert(const Value& value) {
    if (tmp_.capacity()) {
        tmp_.emplace_back(value);
    }
    else {
        tree_->insert(value);
    }
}

template <typename TreeT, typename PayloadT, typename PointT>
void KDTree_eckit<TreeT, PayloadT, PointT>::build() {
    if (tmp_.size()) {
        build(tmp_);
        tmp_.clear();
        tmp_.shrink_to_fit();
    }
}

template <typename TreeT, typename PayloadT, typename PointT>
void KDTree_eckit<TreeT, PayloadT, PointT>::build(std::vector<Value>& values) {
    tree_->build(values);
}


template <typename TreeT, typename PayloadT, typename PointT>
typename KDTree_eckit<TreeT, PayloadT, PointT>::ValueList KDTree_eckit<TreeT, PayloadT, PointT>::do_closestPoints(
    const Point& p, size_t k) const {
    assert_built();
    return tree_->kNearestNeighbours(p, k);
}

template <typename TreeT, typename PayloadT, typename PointT>
typename KDTree_eckit<TreeT, PayloadT, PointT>::Value KDTree_eckit<TreeT, PayloadT, PointT>::do_closestPoint(
    const Point& p) const {
    assert_built();
    return tree_->nearestNeighbour(p);
}

template <typename TreeT, typename PayloadT, typename PointT>
typename KDTree_eckit<TreeT, PayloadT, PointT>::ValueList
KDTree_eckit<TreeT, PayloadT, PointT>::do_closestPointsWithinRadius(const Point& p, double radius) const {
    assert_built();
    return tree_->findInSphere(p, radius);
}

template <typename TreeT, typename PayloadT, typename PointT>
void KDTree_eckit<TreeT, PayloadT, PointT>::assert_built() const {
    if (tmp_.capacity()) {
        throw_AssertionFailed("KDTree was used before calling build()");
    }
}

//------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace util
}  // namespace atlas
