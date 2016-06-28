/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_grids_global_lonlat_ReducedLonLat_h
#define atlas_grids_global_lonlat_ReducedLonLat_h

#include "atlas/grid/global/lonlat/LonLat.h"


namespace atlas {
namespace grid {
namespace global {
namespace lonlat {


/**
 * @brief Reduced LonLat Grid
 *
 * This grid is a special case of the class Structured, with
 * equidistant distribution of latitudes, and a equidistant distribution in zonal
 * direction, which reduce in number going closer towards poles,
 * essentially making the grid more uniform on the sphere
 * It can be constructed with following definition:
 *   N   = number of latitudes in hemisphere
 *   npts_per_lat[] = number of points on each latitude
 */
class ReducedLonLat: public LonLat {

  public:

    static std::string grid_type_str();

    ReducedLonLat(const eckit::Parametrisation&);

    ReducedLonLat(const size_t nlat, const long nlon[], const Domain& dom=Domain::makeGlobal());

    static std::string className();

    virtual eckit::Properties spec() const;

  protected:

    void setup(const eckit::Parametrisation&);

    void setup(const size_t nlat, const long nlon[]);

    void set_typeinfo();

};


} // namespace lonlat
} // namespace global
} // namespace grid
} // namespace atlas


#endif
