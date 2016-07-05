/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_grid_lonlat_RegularLonLat_h
#define atlas_grid_lonlat_RegularLonLat_h

#include "atlas/grid/lonlat/LonLat.h"


namespace atlas {
namespace grid {
namespace lonlat {


class RegularLonLat : public LonLat {

  public:

    static std::string grid_type_str();

    static std::string className();

    RegularLonLat(const eckit::Parametrisation&);

    RegularLonLat(const size_t N, const Domain& dom=Domain::makeGlobal());

    RegularLonLat(const size_t nlon, const size_t nlat, const Domain& dom=Domain::makeGlobal());

    virtual eckit::Properties spec() const;

    size_t nlon() const {
        return Structured::nlon(0);
    }

    double lon(const size_t jlon) const {
        return Structured::lon(0,jlon);
    }

  protected:

    void setup(const eckit::Parametrisation&);

    virtual void set_typeinfo();

};


}  // namespace lonlat
}  // namespace grid
}  // namespace atlas


#endif
