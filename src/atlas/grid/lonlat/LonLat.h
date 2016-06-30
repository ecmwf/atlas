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


class Shift {

  public:
    enum Bits {
        NONE     = 0,
        LON      = (1<<2),
        LAT      = (1<<3)
    };
    Shift() {
        bits_ = NONE;
    }
    Shift(Bits bits) {
        bits_ = static_cast<int>(bits);
    }
    Shift(int  bits) {
        bits_ = bits;
    }
    Shift(bool shift_lon, bool shift_lat) {
        bits_ = NONE;
        if (shift_lon) bits_ |= LON;
        if (shift_lon) bits_ |= LAT;
    }

    bool lat() const      {
        return check(LAT);
    }
    bool lon() const      {
        return check(LON);
    }
    bool lonlat() const   {
        return check(LON|LAT);
    }
    operator bool() const {
        return (bits_ != 0);
    }

  private:
    bool check(int b) const {
        return (bits_ & b) == b;
    }
    int bits_;

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

    bool regular() const {
        return shift_;
    }

  protected:

    void setup(const size_t N, const Domain&);

    void setup(const size_t nlon, const size_t nlat, const Domain&);

    void setup(const long pl[], const size_t nlat, const Domain&);

    virtual void set_typeinfo() = 0;

  protected:

    Shift shift_;

};


}  // namespace lonlat
}  // namespace grid
}  // namespace atlas


#endif
