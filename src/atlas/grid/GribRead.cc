#include "eckit/exception/Exceptions.h"

#include "GribRead.h"

#include "atlas/Mesh.hpp"
#include "atlas/MeshGen.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Parameters.hpp"
#include "atlas/Field.hpp"

#include "Tesselation.h"

//-----------------------------------------------------------------------------

using namespace atlas;

namespace eckit {

//------------------------------------------------------------------------------------------------------

void GribRead::read_nodes_from_grib( grib_handle* h, atlas::Mesh& mesh )
{
    ASSERT( h );

    int err = 0;

    // points to read

    long nb_nodes = 0;
    grib_get_long(h,"numberOfDataPoints",&nb_nodes);

    Tesselation::create_mesh_structure( mesh, nb_nodes );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT(  nodes.bounds()[1] == nb_nodes );

    FieldT<double>& coords  = nodes.field<double>("coordinates");
    FieldT<double>& latlon  = nodes.field<double>("latlon");
    FieldT<int>&    glb_idx = nodes.field<int>("glb_idx");


    if( h == 0 || err != 0 )
        throw std::string("error reading grib");

    grib_iterator *i = grib_iterator_new(h, 0, &err);

    double lat   = 0.;
    double lon   = 0.;
    double value = 0.;

    /// we assume a row first scanning order on the grib
    size_t idx = 0;
    while( grib_iterator_next(i,&lat,&lon,&value) )
    {
        while(lon < 0)    lon += 360;
        while(lon >= 360) lon -= 360;

        glb_idx(idx) = idx;

        latlon(LAT,idx) = lat;
        latlon(LON,idx) = lon;

        atlas::latlon_to_3d( lat, lon, coords.slice(idx) );

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

    atlas::FieldT<double>& field = nodes.create_field<double>(name,1);

    read_field( h, &(field.data()[0]), field.size() );
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

