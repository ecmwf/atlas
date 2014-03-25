#include "eckit/eckit_config.h"
#include "eckit/exception/Exceptions.h"

//------------------------------------------------------------------------------------------------------

#include "atlas/Field.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/Parameters.hpp"

#include "GribWrite.h"

using namespace atlas;

namespace eckit {

//------------------------------------------------------------------------------------------------------

void GribWrite::write( atlas::Mesh& mesh, grib_handle* input_h )
{
    FunctionSpace& nodes     = mesh.function_space( "nodes" );
    FieldT<double>& coords   = nodes.field<double>( "coordinates" );

    const size_t npts = nodes.bounds()[1];

    FieldT<double>& field = nodes.field<double>("field");

    std::vector<double> values (npts);

    ///

    grib_util_grid_spec grid_spec = {0,};

    grid_spec.grid_type = GRIB_UTIL_GRID_SPEC_REGULAR_LL;
    grid_spec.Ni = 60;
    grid_spec.Nj = 30;

    grid_spec.iDirectionIncrementInDegrees = 360. / grid_spec.Ni;
    grid_spec.jDirectionIncrementInDegrees = 180. / (grid_spec.Nj-1);

    grid_spec.longitudeOfFirstGridPointInDegrees =   0. + grid_spec.iDirectionIncrementInDegrees / 2;
    grid_spec.longitudeOfLastGridPointInDegrees  = 360. - grid_spec.iDirectionIncrementInDegrees / 2;

    grid_spec.latitudeOfFirstGridPointInDegrees =   90. - grid_spec.jDirectionIncrementInDegrees / 2;
    grid_spec.latitudeOfLastGridPointInDegrees  =  -90. + grid_spec.jDirectionIncrementInDegrees / 2;

    grib_util_packing_spec pack_spec = {0,};

    pack_spec.packing_type = GRIB_UTIL_PACKING_TYPE_GRID_SIMPLE;
    pack_spec.bitsPerValue = 16;

    int err = 0;

    //    grib_handle* output_h = grib_util_set_spec(input_h,&grid_spec,&pack_spec,0,&(field.data())[0],npts,&err);

    /* flip points for correct grid ordering ???? */

    for ( size_t i = 0; i < grid_spec.Ni; i++ )
    {
        for ( size_t j = 0; j < grid_spec.Nj; j++ )
        {
            size_t idx = (grid_spec.Nj-j-1) + (grid_spec.Ni-i-1)*grid_spec.Nj;
            ASSERT( idx < npts );
            values[idx] = field.data()[j+i*grid_spec.Nj];
        }
    }

    grib_handle* output_h = grib_util_set_spec(input_h,&grid_spec,&pack_spec,0,&(values)[0],npts,&err);

    grib_write_message( output_h, "output.grib", "w" );

    grib_handle_delete(output_h);

#if 0
    const char* filename = "regular_ll_pl_grib1";
    grib_handle* h = grib_handle_new_from_samples(0,"regular_ll_pl_grib1");
    ASSERT( h );

    GRIB_CHECK( grib_set_long(h,"bitsPerValue",16),0 );

    /* set data values*/
    GRIB_CHECK(grib_set_double_array(h,"values",field.data(),npts),0);

    grib_write_message(h,argv[1],"w");

    //    Ni = 16;
    //    Nj = 31;
    GRIB_CHECK( grib_set_long(h,"Ni", Ni ),0 );
    GRIB_CHECK( grib_set_long(h,"Nj", Nj ),0 );

    GRIB_CHECK( grib_set_double(h,"latitudeOfFirstGridPointInDegrees", latitudeOfFirstGridPointInDegrees ),0 );
    GRIB_CHECK( grib_set_double(h,"longitudeOfFirstGridPointInDegrees", longitudeOfFirstGridPointInDegrees ),0 );

    GRIB_CHECK( grib_set_double(h,"latitudeOfLastGridPointInDegrees", latitudeOfLastGridPointInDegrees ),0 );
    GRIB_CHECK( grib_set_double(h,"longitudeOfLastGridPointInDegrees", longitudeOfLastGridPointInDegrees ),0 );

    GRIB_CHECK( grib_set_double(h,"jDirectionIncrementInDegrees", jDirectionIncrementInDegrees ),0 );
    GRIB_CHECK( grib_set_double(h,"iDirectionIncrementInDegrees", iDirectionIncrementInDegrees ),0 );

    free(values);
    grib_handle_delete(h);
#endif
}

//------------------------------------------------------------------------------------------------------

} // namespace eckit

