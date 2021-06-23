/*
 * (C) Crown Copyright 2021, Met Office.
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

#include "atlas/grid/detail/tiles/Tiles.h"

namespace atlas {
namespace cubedspheretiles {

class LFRicCubedSphereTiles : public CubedSphereTiles {
public:

    // constructor
    LFRicCubedSphereTiles( const eckit::Parametrisation& );

    static std::string static_type() { return "LFRicCubedSphereTiles"; }

    virtual std::string type() const override { return static_type(); }

    virtual std::array<std::array<double,6>,2> xy2abOffsets() const override;

    virtual std::array<std::array<double,6>,2> ab2xyOffsets() const override;

    virtual void tile0Rotate( double xyz[] ) const override;

    virtual void tile1Rotate( double xyz[] ) const override;

    virtual void tile2Rotate( double xyz[] ) const override;

    virtual void tile3Rotate( double xyz[] ) const override;

    virtual void tile4Rotate( double xyz[] ) const override;

    virtual void tile5Rotate( double xyz[] ) const override;

    virtual void tile0RotateInverse( double xyz[] ) const override;

    virtual void tile1RotateInverse( double xyz[] ) const override;

    virtual void tile2RotateInverse( double xyz[] ) const override;

    virtual void tile3RotateInverse( double xyz[] ) const override;

    virtual void tile4RotateInverse( double xyz[] ) const override;

    virtual void tile5RotateInverse( double xyz[] ) const override;

    virtual idx_t tileFromXY( const double xy[] ) const override;

    virtual idx_t tileFromLonLat( const double lonlat[] ) const override;

    virtual void enforceXYdomain( double xy[] ) const override;

    virtual void print( std::ostream& ) const override;

private:

};


}  // namespace cubedspheretiles
}  // namespace atlas
