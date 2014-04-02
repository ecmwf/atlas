#ifndef eckit_GribRead_h
#define eckit_GribRead_h

#include "grib_api.h"

#include "atlas/Mesh.hpp"
#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace eckit {

//------------------------------------------------------------------------------------------------------

class GribRead {
public:

    static atlas::grid::Grid* create_grid_from_grib( grib_handle* h );

    static void read_nodes_from_grib( grib_handle* h, atlas::Mesh& mesh );

    static void read_field_from_grib( grib_handle* h, atlas::Mesh& mesh, const std::string& name );

    static void read_field(  grib_handle* h, double* field, size_t size );

};

//------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif

