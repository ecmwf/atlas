/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <string>

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"
#include "atlas/util/Point.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
namespace cubedspheretiles {

class CubedSphereTiles : public util::Object {
public:
    using Spec = util::Config;

public:

    static const CubedSphereTiles* create( );

    static const CubedSphereTiles* create( const eckit::Parametrisation& );

    static const CubedSphereTiles* create( const std::string& );

    virtual std::string type() const = 0;

    virtual std::array<std::array<double,6>,2> xy2abOffsets() const = 0;

    virtual std::array<std::array<double,6>,2> ab2xyOffsets() const = 0;

    virtual void tile0Rotate( double xyz[] ) const = 0;

    virtual void tile1Rotate( double xyz[] ) const = 0;

    virtual void tile2Rotate( double xyz[] ) const = 0;

    virtual void tile3Rotate( double xyz[] ) const = 0;

    virtual void tile4Rotate( double xyz[] ) const = 0;

    virtual void tile5Rotate( double xyz[] ) const = 0;

    virtual void tile0RotateInverse( double xyz[] ) const = 0;

    virtual void tile1RotateInverse( double xyz[] ) const = 0;

    virtual void tile2RotateInverse( double xyz[] ) const = 0;

    virtual void tile3RotateInverse( double xyz[] ) const = 0;

    virtual void tile4RotateInverse( double xyz[] ) const = 0;

    virtual void tile5RotateInverse( double xyz[] ) const = 0;

    virtual idx_t tileFromXY( const double xy[] ) const = 0;

    virtual idx_t tileFromLonLat( const double lonlat[] ) const = 0;

    virtual void enforceXYdomain( double xy[] ) const = 0;

    virtual atlas::PointXY anyXYToFundamentalXY ( const atlas::PointXY & anyXY ) const = 0;

    /// Output to stream
    virtual void print( std::ostream& ) const = 0;

    friend std::ostream& operator<<( std::ostream& s, const CubedSphereTiles& cst ) {
        cst.print( s );
        return s;
    }

};

class FV3CubedSphereTiles;

class LFRicCubedSphereTiles;

} // namespace cubedspheretiles

} // namespace atlas

