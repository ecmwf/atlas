/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file Polygon.h
/// @author Pedro Maciel
/// @author Willem Deconinck
/// @date September 2017

#pragma once

#include <iosfwd>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "atlas/util/mdspan.h"

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"
#include "atlas/util/Point.h"
#include "atlas/util/VectorOfAbstract.h"

namespace eckit {
class PathName;
}

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

/// Polygon
class Polygon : public std::vector<idx_t> {
public:
    // -- Types

    struct edge_t : std::pair<idx_t, idx_t> {
        edge_t(idx_t A, idx_t B): std::pair<idx_t, idx_t>(A, B) {}

        edge_t reverse() const { return edge_t(std::pair<idx_t, idx_t>::second, std::pair<idx_t, idx_t>::first); }

        struct LessThan {
            bool operator()(const edge_t& e1, const edge_t& e2) const {
                // order ascending by 'first'
                return (e1.first < e2.first ? true : e1.first > e2.first ? false : e1.second < e2.second);
            }
        };
    };

    using edge_set_t  = std::set<edge_t, typename edge_t::LessThan>;
    using container_t = std::vector<idx_t>;

    // -- Constructors

    Polygon();
    Polygon(const edge_set_t&);

    // -- Operators

    operator bool() const;

    Polygon& operator+=(const Polygon&);

    // -- Methods

    void print(std::ostream&) const;

    // -- Friends

    friend std::ostream& operator<<(std::ostream& s, const Polygon& p) {
        p.print(s);
        return s;
    }

protected:
    void setup(const edge_set_t&);
};

//------------------------------------------------------------------------------------------------------

class PartitionPolygons;

class PartitionPolygon : public Polygon, util::Object {
public:
    using Polygon::Polygon;
    using PointsXY     = std::vector<Point2>;
    using PointsLonLat = std::vector<Point2>;
    class RectangleXY {
    public:
        RectangleXY() = default;
        RectangleXY(const std::array<double,2>& x, const std::array<double,2>& y) :
            xmin_{x[0]}, xmax_{x[1]}, ymin_{y[0]}, ymax_{y[1]} {}

        double xmin() const { return xmin_; }
        double xmax() const { return xmax_; }
        double ymin() const { return ymin_; }
        double ymax() const { return ymax_; }
    private:
        double xmin_{std::numeric_limits<double>::max()};
        double xmax_{std::numeric_limits<double>::max()};
        double ymin_{std::numeric_limits<double>::max()};
        double ymax_{std::numeric_limits<double>::max()};
    };

    /// @brief Return inscribed rectangular domain (not rotated)
    virtual const RectangleXY& inscribedDomain() const;

    /// @brief Return value of halo
    virtual idx_t halo() const { return 0; }

    /// @brief Return the memory footprint of the Polygon
    virtual size_t footprint() const { return 0; }

    /// @brief Output a python script that plots the partition
    virtual void outputPythonScript(const eckit::PathName&, const eckit::Configuration& = util::NoConfig()) const {}

    /// @brief Output a JSON file with partition polygons
    virtual std::string json(const eckit::Configuration& = util::NoConfig()) const;

    /// @brief All (x,y) coordinates defining a polygon. Last point should match first.
    virtual PointsXY xy() const = 0;

    /// @brief All (lon,lat) coordinates defining a polygon. Last point should match first.
    virtual PointsLonLat lonlat() const = 0;

    virtual void allGather(PartitionPolygons&) const = 0;
};

//------------------------------------------------------------------------------------------------------

class ExplicitPartitionPolygon : public util::PartitionPolygon {
public:
    explicit ExplicitPartitionPolygon(PointsXY&& points):
        ExplicitPartitionPolygon(std::move(points), RectangleXY()) {}

    explicit ExplicitPartitionPolygon(PointsXY&& points, const RectangleXY& inscribed):
        points_(std::move(points)), inscribed_(inscribed) {
        setup(compute_edges(points_.size()));
    }

    PointsXY xy() const override { return points_; }
    PointsLonLat lonlat() const override { return points_; }

    void allGather(util::PartitionPolygons&) const override;

    const RectangleXY& inscribedDomain() const override { return inscribed_; }

private:
    static util::Polygon::edge_set_t compute_edges(idx_t points_size);


private:
    std::vector<Point2> points_;
    RectangleXY inscribed_;
};  // namespace util

//------------------------------------------------------------------------------------------------------

class PartitionPolygons : public VectorOfAbstract<PartitionPolygon> {
public:
    std::string json(const eckit::Configuration& = util::NoConfig()) const;
};

//------------------------------------------------------------------------------------------------------

class PolygonCoordinates {
public:
    using Vector = VectorOfAbstract<PolygonCoordinates>;
    // -- Constructors

    PolygonCoordinates(const Polygon&, const double x[], const double y[], size_t xstride, size_t ystride, bool removeAlignedPoints);

    template <typename Extents, typename LayoutPolicy, typename AccessorPolicy>
    PolygonCoordinates(const Polygon& poly, atlas::mdspan<const double, Extents, LayoutPolicy, AccessorPolicy> coordinates, bool removeAlignedPoints): 
        PolygonCoordinates(poly, &coordinates(0,0), &coordinates(0,1), coordinates.stride(1), coordinates.stride(1), removeAlignedPoints ) {}

    template <typename PointContainer>
    PolygonCoordinates(const PointContainer& points);

    template <typename PointContainer>
    PolygonCoordinates(const PointContainer& points, bool removeAlignedPoints);

    // -- Destructor

    virtual ~PolygonCoordinates();

    // -- Methods

    /// @brief Point-in-partition test
    /// @param[in] P given point
    /// @return if point is in polygon
    virtual bool contains(const Point2& P) const = 0;

    const Point2& coordinatesMax() const;
    const Point2& coordinatesMin() const;
    const Point2& centroid() const;

    template <typename Index>
    const Point2& operator[](Index i) const {
        return coordinates_[i];
    }

    idx_t size() const { return coordinates_.size(); }

    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& out, const PolygonCoordinates& pc);

protected:
    // -- Members

    Point2 coordinatesMin_;
    Point2 coordinatesMax_;
    Point2 centroid_;
    std::vector<Point2> coordinates_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
