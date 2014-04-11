#include "eckit/exception/Exceptions.h"
#include "eckit/log/Timer.h"
#include "eckit/log/Seconds.h"
#include "eckit/geometry/Point3.h"
#include "eckit/grib/GribAccessor.h"

#include "atlas/Mesh.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Parameters.hpp"
#include "atlas/Field.hpp"

#include "atlas/grid/GribRead.h"
#include "atlas/grid/Unstructured.h"
#include "atlas/grid/Tesselation.h"
#include "atlas/grid/RegularLatLonGrid.h"

//-----------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::grid;

namespace eckit { /// @todo this is still in eckit namespace because we plan to move it back to eckit::grib

//------------------------------------------------------------------------------------------------------

typedef std::vector< grid::Grid::Point > PointList;

static PointList* read_number_of_data_points(grib_handle *h)
{
   // points to read
   long nb_nodes = 0;
   grib_get_long(h,"numberOfDataPoints",&nb_nodes);

   int err = 0;
   grib_iterator *i = grib_iterator_new(h, 0, &err);
   if(err != 0 )
      throw std::string("error reading grib");

   PointList* pts = new PointList(nb_nodes);

   double lat   = 0.;
   double lon   = 0.;
   double value = 0.;

   size_t idx = 0;
   while( grib_iterator_next(i,&lat,&lon,&value) )
   {
      (*pts)[idx].assign(lat,lon);
      ++idx;
   }
   grib_iterator_delete(i);

   ASSERT( idx == nb_nodes );
   return pts;
}

grid::Grid* GribRead::create_grid_from_grib(grib_handle *h)
{
   ASSERT( h );
   if ( 0 == h) throw std::string("GribRead::create_grid_from_grib NULL grib_handle");

   char string_value[64];
   size_t len = sizeof(string_value)/sizeof(char);
   if (grib_get_string(h,"gridType",string_value,&len)  != 0) {
      throw std::runtime_error("grib_get_string failed for gridType") ;
   }

   if (strncasecmp(string_value,"regular_ll",10) == 0) {
      // Custom arguments for Lat long grid
      return new grid::RegularLatLonGrid( h );
   }
//   else if (strncasecmp(string_value,"sh",2) == 0) {
//      return new grid::SphericalHarmonicGrid( pts, grib_hash(h) );
//   }
//   else if (strncasecmp(string_value,"reduced_ll",10) == 0) {
//      return new grid::ReducedLatLonGrid( pts, grib_hash(h) );
//   }
//   else if (strncasecmp(string_value,"reduced_gg",10) == 0) {
//      return new grid::ReducedGuassianGrid( pts, grib_hash(h) );
//   }
//   else if (strncasecmp(string_value,"regular_gg",10) == 0) {
//      return new grid::RegularGuassianGrid( pts, grib_hash(h) );
//   }

   // Unknown grid type, get extract data points form the grib handle
   return new grid::Unstructured( read_number_of_data_points(h), grib_hash(h) );
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

    ASSERT(  nodes.bounds()[1] == nb_nodes );

    FieldT<double>& coords  = nodes.field<double>("coordinates");
    FieldT<double>& latlon  = nodes.field<double>("latlon");
    FieldT<int>&    glb_idx = nodes.field<int>("glb_idx");

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

        latlon(LAT,idx) = lat;
        latlon(LON,idx) = lon;

        eckit::geometry::latlon_to_3d( lat, lon, coords.slice(idx) );

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

    read_field( h, &field.data()[0], field.size() );
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


void GribRead::known_grid_types(std::set<std::string>& grids)
{
   grids.insert("regular_ll");
   grids.insert("reduced_ll");
   grids.insert("mercator");
   grids.insert("lambert");
   grids.insert("polar_stereographic");
   grids.insert("UTM");
   grids.insert("simple_polyconic");
   grids.insert("albers");
   grids.insert("miller");
   grids.insert("rotated_ll");
   grids.insert("stretched_ll");
   grids.insert("stretched_rotated_ll");
   grids.insert("regular_gg");
   grids.insert("rotated_gg");
   grids.insert("stretched_gg");
   grids.insert("stretched_rotated_gg");
   grids.insert("reduced_gg");
   grids.insert("sh");
   grids.insert("rotated_sh");
   grids.insert("stretched_sh");
   grids.insert("stretched_rotated_sh");
   grids.insert("space_view");
   grids.insert("unknown");
   grids.insert("unknown_PLPresent");
}

//---------------------------------------------------------------------------------------------------------

} // namespace eckit

