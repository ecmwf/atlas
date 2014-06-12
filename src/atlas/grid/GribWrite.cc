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

#include "eckit/eckit_config.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/grib/GribHandle.h"
#include "eckit/utils/Translator.h"
#include "eckit/memory/ScopedPtr.h"
#include "eckit/io/DataHandle.h"

#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/Parameters.hpp"

#include "atlas/grid/FieldSet.h"
#include "atlas/grid/GribWrite.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using namespace atlas::grid;

#define DBG     std::cout << Here() << std::endl;
#define DBGX(x) std::cout << #x << " -> " << x << std::endl;

namespace atlas {

//------------------------------------------------------------------------------------------------------

GribHandle *GribWrite::create_handle(const Grid &)
{
    NOTIMP;
}

void GribWrite::write( const FieldSet& fields, const PathName& opath )
{
    for( size_t i = 0; i < fields.size(); ++i )
    {
        PathName pi( opath.asString() + "." + Translator<size_t,std::string>()(i) );
        GribWrite::write(fields[i], pi);
    }
}

void GribWrite::clone( const FieldSet& fields, const PathName& src, const PathName& opath )
{
    bool overwrite = true;

    if( opath.exists() )
        opath.unlink();

    eckit::ScopedPtr<DataHandle> of( opath.fileHandle(overwrite) ); AutoClose of_close(*of);

    ASSERT(of);

    of->openForWrite(0);

    for( size_t i = 0; i < fields.size(); ++i )
    {
        GribWrite::clone(fields[i], src, *of);
    }
}

void GribWrite::write(const FieldHandle &field, const PathName &opath)
{
    NOTIMP;
}

void GribWrite::clone(const FieldHandle& field, const PathName& gridsec, DataHandle& out )
{
    FILE* fh = ::fopen( gridsec.asString().c_str(), "r" );
    if( fh == 0 )
        throw ReadError( std::string("error opening file ") + gridsec );

    int err = 0;
    grib_handle* clone_h = grib_handle_new_from_file(0,fh,&err);
    if( clone_h == 0 || err != 0 )
        throw ReadError( std::string("error reading grib file ") + gridsec );

    GribHandle ch(clone_h);
    ScopedPtr<GribHandle> h( GribWrite::clone( field, ch ) );

    //    GRIB_CHECK( grib_write_message(h->raw(),fname.asString().c_str(),"w"), 0 );

    // dump the handle to the DataHandle
    const void* buffer = NULL;
    size_t size = 0;

    GRIB_CHECK( grib_get_message( h->raw(), &buffer, &size), 0);

    out.write(buffer, size);
}

GribHandle* GribWrite::clone(const FieldHandle& field, GribHandle& gridsec )
{
    const Field& f = field.data();
    const size_t npts = f.size();

    // check number of points matches

    long nb_nodes = 0;
    GRIB_CHECK( grib_get_long(gridsec.raw(),"numberOfDataPoints",&nb_nodes), 0 );
    ASSERT( npts == f.size() );

    GribHandle& meta = field.grib();

    ///@todo move this to the eckit::grib interface
    int err=0;
    int what = GRIB_SECTION_GRID;
    grib_handle* h = grib_util_sections_copy( gridsec.raw(), meta.raw(), what, &err); GRIB_CHECK(err,"grib_util_sections_copy()");

    ///@todo move this to the eckit::grib interface
    GRIB_CHECK( grib_set_double_array(h, "values", f.data<double>(),npts), 0 );

    return new GribHandle(h);
}

#if 0
void GribWrite::write( atlas::grid::FieldHandle& field, grib_handle* input_h )
{
    FieldT<double>& f = field.data();
    const size_t npts = f.size();

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
            values[idx] = f[j+i*grid_spec.Nj];
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

    grib_handle_delete(h);
#endif
}
#endif

//------------------------------------------------------------------------------------------------------

} // namespace atlas

