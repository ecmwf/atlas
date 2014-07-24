/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grid_GridSpec_H
#define atlas_grid_GridSpec_H

#include "eckit/value/Properties.h"

#include "atlas/grid/Grid.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------


/// A concrete class that holds specification that uniquely identifies a Grid
/// The description of the grid is added as name value pairs
/// This class will provides a short name for a GRID (i.e QG48_1)
/// This allows for easier matching with samples files.
/// However this interface is independent of GRIB/NETCDF, because:
///
///      DECODE                   ATLAS                      ENCODE
///  NetCDFParams ---->|-------|         |----------|------>NetCDFGridWrite
///                    | Grid  |<------> | GridSpec |
///  GribParams ------>|-------|         |----------|------>GribGridWrite
///
/// Uses default copy constructor, assignment and equality operators


class GridSpec : public eckit::Properties {
public:

    GridSpec(const std::string& grid_type);

    ~GridSpec();

    /// returns the gridType. currently this matches grid _type found in GRIB
    std::string grid_type() const;

    void uid(const std::string&);
    std::string uid() const;

    /// Helper functions, Used to build up a grid spec
    void set_latitudes(const std::vector<double>& latitudes);
    void set_rgspec(const std::vector<long>&  rgSpec);
    void set_bounding_box(const Grid::BoundBox& bbox );

    void get_latitudes(std::vector<double>& latitudes) const;
    void get_rgspec(std::vector<long>& rgSpec) const;
    void get_bounding_box(Grid::BoundBox& bbox ) const;

    friend std::ostream& operator<<( std::ostream& os, const GridSpec& v) { v.print(os); return os;}

private:

    void print( std::ostream& ) const;

};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
