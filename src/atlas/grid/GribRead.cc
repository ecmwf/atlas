/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"
#include "eckit/log/Timer.h"
#include "eckit/log/Seconds.h"
#include "eckit/geometry/Point3.h"
#include "eckit/grib/GribAccessor.h"
#include "eckit/grib/GribHandle.h"

#include "atlas/mesh/Mesh.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/mesh/Field.hpp"

#include "atlas/grid/GribRead.h"
#include "atlas/grid/Tesselation.h"
#include "atlas/grid/GribGridBuilder.h"

//-----------------------------------------------------------------------------

using namespace std;
using namespace atlas;
using namespace atlas::grid;

namespace eckit { /// @todo this is still in eckit namespace because we plan to move it back to eckit::grib

//------------------------------------------------------------------------------------------------------

grid::Grid::Ptr GribRead::create_grid_from_grib( grib_handle* h )
{
   if ( !h )
       throw std::runtime_error("GribRead::create_grid_from_grib NULL grib_handle");

   return GRIBGridBuilder::instance().build_grid_from_grib_handle( h );
}

void GribRead::read_nodes_from_grib( grib_handle* h, atlas::Mesh& mesh )
{
    ASSERT( h );

    int err = 0;

    // points to read

    long nb_nodes = 0;
    grib_get_long(h,"numberOfDataPoints",&nb_nodes);

    Tesselation::create_mesh_structure( mesh, nb_nodes );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT(  nodes.extents()[0] == nb_nodes );

    ArrayView<double,2> coords  ( nodes.field("coordinates") );
    ArrayView<double,2> latlon  ( nodes.field("latlon") );
    ArrayView<int,   1> glb_idx ( nodes.field("glb_idx") );

    grib_iterator *i = grib_iterator_new(h, 0, &err);

    if( h == 0 || err != 0 )
        throw std::string("error reading grib");

    double lat   = 0.;
    double lon   = 0.;
    double value = 0.;

//    Timer t("inside read_nodes_from_grib");

    /// we assume a row first scanning order on the grib
    size_t idx = 0;
    while( grib_iterator_next(i,&lat,&lon,&value) )
    {
        while(lon < 0)    lon += 360;
        while(lon >= 360) lon -= 360;

        glb_idx(idx) = idx;

        latlon(idx,LAT) = lat;
        latlon(idx,LON) = lon;

        eckit::geometry::latlon_to_3d( lat, lon, coords[idx].data() );

        ++idx;
    }
    grib_iterator_delete(i);

    ASSERT( idx == nb_nodes );
}

//------------------------------------------------------------------------------------------------------

void GribRead::read_field_from_grib(  grib_handle* h, atlas::Mesh& mesh, const std::string& name )
{
    ASSERT( h );
    ASSERT( mesh.has_function_space("nodes") );

    atlas::FunctionSpace& nodes  = mesh.function_space( "nodes" );

    atlas::Field& field = nodes.create_field<double>(name,1);

    read_field( h, field.data<double>(), field.size() );
}

//------------------------------------------------------------------------------------------------------

void GribRead::read_field(  grib_handle* h, double* field, size_t size )
{    
    ASSERT( h );

    long nb_nodes = 0;
    grib_get_long(h,"numberOfDataPoints",&nb_nodes);

    if( nb_nodes != size )
    {
        std::ostringstream msg;
        msg << "number of data points in grib " << nb_nodes
            << " differs from field " << size
            << std::endl;
        throw SeriousBug( msg.str() );
    }

    int err = 0;
    grib_iterator *i = grib_iterator_new(h, 0, &err);

    double lat   = 0.;
    double lon   = 0.;
    double value = 0.;

    size_t in = 0;
    while(grib_iterator_next(i,&lat,&lon,&value))
    {
        if( in >= size )
            throw SeriousBug( "field is of incorrect size -- too many points" );
        field[in] = value;
        ++in;
    }
    grib_iterator_delete(i);

    if( in != size )
        throw SeriousBug( "field is of incorrect size -- too little points" );

}


//---------------------------------------------------------------------------------------------------------

} // namespace eckit

