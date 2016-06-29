/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_grids_global_lonlat_ShiftedLat_h
#define atlas_grids_global_lonlat_ShiftedLat_h

#include "atlas/grid/global/lonlat/LonLat.h"


namespace atlas {
namespace grid {
namespace global {
namespace lonlat {


class ShiftedLat: public LonLat {

  public:

    static std::string grid_type_str();

    static std::string className();

    ShiftedLat(const eckit::Parametrisation&);

    ShiftedLat(const size_t N);

    ShiftedLat(const size_t nlon, const size_t nlat);

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
}  // namespace global
}  // namespace grid
}  // namespace atlas


#endif
