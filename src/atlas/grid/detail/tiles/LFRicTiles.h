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

    virtual idx_t tileFromXY( const double xy[] ) const override;

    virtual idx_t tileFromLonLat( const double lonlat[] ) const override;

    virtual void enforceXYdomain( double xy[] ) const override;

    virtual void print( std::ostream& ) const override;


private:


};


}  // namespace cubedspheretiles
}  // namespace atlas
