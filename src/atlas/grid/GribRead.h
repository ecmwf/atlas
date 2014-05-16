#ifndef eckit_GribRead_h
#define eckit_GribRead_h
/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "grib_api.h"
#include <vector>

#include "atlas/mesh/Mesh.hpp"
#include "atlas/grid/Grid.h"
#include "eckit/geometry/Point3.h"

//-----------------------------------------------------------------------------

namespace eckit {

//------------------------------------------------------------------------------------------------------

class GribRead {
public:

    static atlas::grid::Grid* create_grid_from_grib( grib_handle* h );

    static void read_nodes_from_grib( grib_handle* h, atlas::Mesh& mesh );

    static void read_field_from_grib( grib_handle* h, atlas::Mesh& mesh, const std::string& name );

    static void read_field(  grib_handle* h, double* field, size_t size );

    /// Returns the list of all known grib grid types.
    /// Allow better error handling during Grid construction.
    /// This really belongs in grib_api,
    static void known_grid_types(std::set<std::string>&);

    // The return list is new allocated,
    typedef std::vector< atlas::grid::Grid::Point > PointList;
    static PointList* read_number_of_data_points(grib_handle *h);

    static void comparePointList(const std::vector<atlas::grid::Grid::Point>&, double epsilon,grib_handle*);

};

//------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif

