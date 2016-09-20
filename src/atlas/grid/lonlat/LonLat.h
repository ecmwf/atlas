/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_grid_lonlat_LonLat_h
#define atlas_grid_lonlat_LonLat_h

#include "atlas/grid/Structured.h"


namespace atlas {
namespace grid {
namespace lonlat {


struct Shift {

    enum Bits {
        NONE = 0,
        LAT  = (1<<1),
        LON  = (1<<2)
    };

    Shift(int bits=NONE) : bits_(bits) {
    }

    Shift(bool shift_lon, bool shift_lat) : bits_((shift_lon? LON:NONE) | (shift_lat? LAT:NONE)) {
    }

    bool operator()(int bits) const {
        return (bits_ & bits) == bits;
    }

    const int bits_;

};


/**
 * @brief (Structured) LonLat Grid
 *
 * This grid is a special case of the class Structured, with
 * equidistant distribution of latitudes, and a equidistant distribution in zonal
 * direction, which reduce in number going closer towards poles,
 * essentially making the grid more uniform on the sphere
 * It can be constructed with following definition:
 * * N   = number of latitudes in hemisphere
 * * npts_per_lat[] = number of points on each latitude 
 */
class LonLat: public Structured {

  public:

    static std::string grid_type_str();

    static std::string className();

    LonLat(const Shift&, const Domain& dom=Domain::makeGlobal());

    const Shift& shifted() const {
        return shift_;
    }

    virtual const Domain& domain() const {
        return domain_;
    }
    
    //virtual const Projection& projection() const {
    //		return projection_;
    //}


  protected:

    void setup(const size_t N, const Domain&);

    void setup(const size_t nlon, const size_t nlat, const Domain&);

    void setup(const long pl[], const size_t nlat, const Domain&);

    virtual void set_typeinfo() = 0;

    static eckit::Value domain_spec(const Domain& dom);

  protected:

    Shift shift_;

    Domain domain_;
    
    //Projection & projection_;

};


}  // namespace lonlat
}  // namespace grid
}  // namespace atlas


#endif
