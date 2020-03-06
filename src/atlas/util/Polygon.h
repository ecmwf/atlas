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
#include <set>
#include <utility>
#include <vector>

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"
#include "atlas/util/Point.h"

namespace eckit {
class PathName;
}

namespace atlas {
class Field;
class RectangularDomain;
}  // namespace atlas

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

/// Polygon
class Polygon : public std::vector<idx_t> {
public:
    // -- Types

    struct edge_t : std::pair<idx_t, idx_t> {
        edge_t( idx_t A, idx_t B ) : std::pair<idx_t, idx_t>( A, B ) {}

        edge_t reverse() const { return edge_t( std::pair<idx_t, idx_t>::second, std::pair<idx_t, idx_t>::first ); }

        struct LessThan {
            bool operator()( const edge_t& e1, const edge_t& e2 ) const {
                // order ascending by 'first'
                return ( e1.first < e2.first ? true : e1.first > e2.first ? false : e1.second < e2.second );
            }
        };
    };

    using edge_set_t  = std::set<edge_t, typename edge_t::LessThan>;
    using container_t = std::vector<idx_t>;

    // -- Constructors

    Polygon();
    Polygon( const edge_set_t& );

    // -- Operators

    operator bool() const;

    Polygon& operator+=( const Polygon& );

    // -- Methods

    void print( std::ostream& ) const;

    // -- Friends

    friend std::ostream& operator<<( std::ostream& s, const Polygon& p ) {
        p.print( s );
        return s;
    }

protected:
    void setup( const edge_set_t& );
};

//------------------------------------------------------------------------------------------------------

class PartitionPolygon : public Polygon, util::Object {
public:
    using Polygon::Polygon;

    /// @brief Return inscribed rectangular domain (not rotated)
    virtual const RectangularDomain& inscribedDomain() const;

    /// @brief Return the memory footprint of the Polygon
    virtual idx_t halo() const { return 0; }

    /// @brief Return the memory footprint of the Polygon
    virtual size_t footprint() const { return 0; }

    /// @brief Output a python script that plots the partition
    virtual void outputPythonScript( const eckit::PathName&, const eckit::Configuration& = util::NoConfig() ) const {}

    virtual const std::vector<Point2>& xy() const = 0;

    virtual const std::vector<Point2>& lonlat() const = 0;
};

//------------------------------------------------------------------------------------------------------

class PolygonCoordinates {
public:
    // -- Constructors

    PolygonCoordinates( const Polygon&, const atlas::Field& coordinates, bool removeAlignedPoints );

    template <typename PointContainer>
    PolygonCoordinates( const PointContainer& points );

    template <typename PointContainer>
    PolygonCoordinates( const PointContainer& points, bool removeAlignedPoints );

    // -- Destructor

    virtual ~PolygonCoordinates();

    // -- Methods

    /*
   * Point-in-partition test
   * @param[in] P given point
   * @return if point is in polygon
   */
    virtual bool contains( const Point2& P ) const = 0;

    const Point2& coordinatesMax() const;
    const Point2& coordinatesMin() const;

protected:
    // -- Members

    Point2 coordinatesMin_;
    Point2 coordinatesMax_;
    std::vector<Point2> coordinates_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
