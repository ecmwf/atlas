/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/atlas_config.h"

#ifdef ECKIT_HAVE_GRIB

#include "grib_api.h" // remove this once we use only eckit::grib

#include "eckit/eckit_config.h"
#include "eckit/config/Resource.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/utils/Translator.h"
#include "eckit/memory/ScopedPtr.h"
#include "eckit/io/DataHandle.h"
#include "eckit/io/FileHandle.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/parser/StringTools.h"

#include "eckit/grib/GribParams.h"
#include "eckit/grib/GribHandle.h"
#include "eckit/grib/GribMutator.h"

#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Parameters.h"

#include "atlas/FieldSet.h"
#include "atlas/io/Grib.h"
#include "atlas/GridSpec.h"

//------------------------------------------------------------------------------------------------------

using namespace std;
using namespace eckit;
using namespace eckit::grib;
using namespace atlas;

namespace atlas {
namespace io {

static std::string map_uid_to_grib_sample_file(const std::string& short_name, long edition);

//------------------------------------------------------------------------------------------------------

Grid::Ptr Grib::create_grid(GribHandle& gh)
{
	GribParams* gp = GribParams::create(gh);

	ASSERT( gp );

	return Grid::create( *gp );
}

GribHandle::Ptr Grib::create_handle( const Grid& grid, long edition )
{
    // determine choice of editionNumber from a resource

    if( edition == 0 )
        edition = Resource<unsigned>( "NewGribEditionNumber", 2 );

    // From the Grid get the Grid Spec
	GridSpec grid_spec = grid.spec();

    grib_handle* gh = 0;
    std::string sample_file;

    // first match GridSpec uid, directly to a samples file

	sample_file = map_uid_to_grib_sample_file( grid_spec.uid(), edition );
    if( !sample_file.empty() )
    {
        gh = grib_handle_new_from_samples(0,sample_file.c_str() );
        if( !gh )
            throw SeriousBug( "Failed to create GribHandle from sample: " + sample_file, Here() );
    }
    else // if this fails, then try looking on disk
    {
		sample_file = Grib::grib_sample_file(grid_spec,edition);
        if (!sample_file.empty())
        {
            gh = grib_handle_new_from_samples(0,sample_file.c_str() );
            if( !gh )
                throw SeriousBug( "Failed to create GribHandle from sample: " + sample_file, Here() );
        }
    }

    if(!gh)
        throw SeriousBug( "Failed to create GribHandle from Grib", Here() );

	GribHandle::Ptr gh_ptr( new GribHandle(gh) );



	write_gridspec_to_grib( grid.spec(), *gh_ptr );

	return gh_ptr;
}


void Grib::determine_grib_samples_dir(std::vector<std::string>& sample_paths)
{
   // The grib samples has the format: <dir-path1> : <dir-path2> : <dir-path3>
   // Also only grib files in these directories that have a '.tmpl' are considered
   //
   char* paths = NULL;
#ifdef HAVE_GRIB_API_1130
   paths = grib_samples_path(NULL);
#endif

   if (paths)
   {
      // Expect <path1>:<path2>:<path3:
      sample_paths = StringTools::split(":", std::string(paths));
      return;
   }

   // GRIB_SAMPLES_PATH may be set by the user, if set use in preference to GRIB_API_PATH
   char* samples_dir = getenv("GRIB_SAMPLES_PATH");
   if (!samples_dir) {

      // GRIB_API_PATH is set for all grib versions at ecmwf, when using modules functionality
      samples_dir = getenv("GRIB_API_PATH");
      if (samples_dir) {
         std::string path = samples_dir;
         path += "/share/samples";
         sample_paths.push_back(path);
         return;
      }
	   throw SeriousBug(string("GRIB_SAMPLES_PATH not defined"),Here()) ;
   }

   sample_paths = StringTools::split(":", std::string(samples_dir));
}

bool match_grid_spec_with_sample_file( const GridSpec& g_spec,
									   GribHandle& handle,
									   long edition,
									   const std::string& file_path )
{
	if ( g_spec.grid_type() != handle.gridType() ) {
        //Log::info() << "grid_type in GridSpec " << g_spec.grid_type() << " does not match " << grib_grid_type << " in samples file " << file_path << " IGNORING " << std::endl;
        return false;
    }

//    eckit::Value spec_nj = g_spec.get("Nj");
//    if (!spec_nj.isNil()) {
//        long spec_nj = spec_nj;
//        long grib_nj = 0;
//        if (grib_get_long(handle.raw(),"Nj",&grib_nj) == 0 ) {
//            if (spec_nj != grib_nj ) {
//                //Log::info() << "Grib::match_grid_spec_with_sample_file, Nj in GridSpec " << spec_nj << " does not match  " << grib_nj << " in samples file " << file_path << " IGNORING " << std::endl;
//                return false;
//            }
//        }
//    }
//    eckit::Value spec_ni = g_spec.get("Ni");
//    if (!spec_ni.isNil()) {
//        long spec_ni = spec_ni;
//        long grib_ni = 0;
//        if (grib_get_long(handle.raw(),"Ni",&grib_ni) == 0 ) {
//            if (spec_ni != grib_ni ) {
//                //Log::info() << "Grib::match_grid_spec_with_sample_file, Ni in GridSpec " << spec_ni << " does not match  " << grib_ni << " in samples file " << file_path << " IGNORING " << std::endl;
//                return false;
//            }
//        }
//    }

	if( handle.edition() != edition ) {
		//Log::info() << "Grib::match_grid_spec_with_sample_file, edition_number passed in " << edition << " does not match grib" << edition << " in samples file " << file_path << " IGNORING " << std::endl;
        return false;
    }

    return true;
}

std::string Grib::grib_sample_file( const GridSpec& g_spec, long edition )
{
    // Note: many of the grib samples files are not UNIQUE in their grid specification:
    // i.e
    //   GRIB2.tmpl                        -> GridSpec[ regular_ll, LL31_16_2, Ni:16, Nj:31, typeOfLevel:surface ]
    //   regular_ll_pl_grib2.tmpl          -> GridSpec[ regular_ll, LL31_16_2, Ni:16, Nj:31 ]
    //   regular_ll_sfc_grib2.tmpl         -> GridSpec[ regular_ll, LL31_16_2, Ni:16, Nj:31 ]
    //
    //   reduced_gg_ml_grib1               -> GridSpec[ reduced_gg, QG32_1, Nj:64 ]
    //   reduced_gg_pl_32_grib1            -> GridSpec[ reduced_gg, QG32_1, Nj:64 ]
    //   reduced_gg_ml_grib2               -> GridSpec[ reduced_gg, QG32_2, Nj:64 ]
    //   reduced_gg_pl_32_grib2            -> GridSpec[ reduced_gg, QG32_2, Nj:64 ]
    //
    // Others are just plain wrong, i.e
    //   polar_stereographic_pl_grib2.tmpl -> GridSpec[ rotated_ll, RL31_2, Ni:16, Nj:31, editionNumber:2 ]
    //
    // From the grid spec, we will look at the grid samples, and find the closest match

    std::vector<std::string> sample_paths;
    determine_grib_samples_dir(sample_paths);

    if ( sample_paths.empty() )
        throw SeriousBug(string("Error no sample paths found"),Here()) ;

    for(size_t path = 0; path < sample_paths.size(); ++path)
    {
        std::string grib_samples_dir = sample_paths[path];

        if (grib_samples_dir.empty())
            throw SeriousBug(string("Error, empty samples path. Could not create handle from grid"),Here()) ;

        PathName dir_path(grib_samples_dir);

        if( !dir_path.exists() ) continue;
        if( !dir_path.isDir()  ) continue;

        std::vector<PathName> files;
        std::vector<PathName> directories;
        dir_path.children(files,directories);

        for(size_t i = 0; i < files.size(); i++)
        {
           try
           {
              GribHandle grib_h( files[i].localPath() );

              std::string fname = files[i].localPath();

              if( match_grid_spec_with_sample_file(g_spec,grib_h,edition,fname))
              {
                 // remove .tmpl extension
                 eckit::LocalPathName path(fname);
                 LocalPathName base_name = path.baseName(false);
                 string grib_sample_file = base_name.localPath();
                 return grib_sample_file;
              }
           }
           catch ( const std::exception & ex )
           {
              Log::info() << files[i].localPath() << " " << ex.what() << std::endl;
           }
        }
    }

    Log::info() << "Could find grib samples match for grid_spec " << g_spec << std::endl;
    return std::string();
}

static std::string map_uid_to_grib_sample_file(const std::string& uid, long edition)
{
    using std::string;

    long ns[14] = {32,48,80,128,160,200,256,320,400,512,640,1024,1280,2000};

    std::map<std::string,std::string> uid_to_sample;

    for( size_t i = 0; i < sizeof(ns)/sizeof(long); ++i)
        uid_to_sample[ "reduced_gg_" + Translator<long,string>()(ns[i]) ] = string("reduced_gg_pl_" + Translator<long,string>()(ns[i]) );

    string r;

    std::map<string,string>::const_iterator i = uid_to_sample.find(uid);
    if (i != uid_to_sample.end())
    {
      r = (*i).second + "_grib" + Translator<long,string>()(edition);
    }

    return r;
}

void Grib::write( const FieldSet& fields, const PathName& opath )
{
    for( size_t i = 0; i < fields.size(); ++i )
    {
        PathName pi( opath.asString() + "." + Translator<size_t,std::string>()(i) );
		Grib::write(fields[i], pi);
    }
}

void Grib::write(const Field& fh, DataHandle& dh)
{
	GribHandle::Ptr gh = Grib::create_handle( fh.grid(), fh.grib().edition() );

    if( !gh )
		throw SeriousBug("Failed to create GribHandle from Field", Here());

    if( !gh->raw() )
		throw SeriousBug("Failed to create GribHandle from Field", Here());

    GribHandle::Ptr h = clone(fh,*gh);

    // dump the handle to the DataHandle
    const void* buffer = NULL;
    size_t size = 0;

    GRIB_CHECK( grib_get_message( h->raw(), &buffer, &size), 0);

    dh.write(buffer, size);
}

GribHandle::Ptr Grib::write(const Field& fh)
{
	GribHandle::Ptr gh = Grib::create_handle( fh.grid(), fh.grib().edition() );

    if( !gh )
		throw SeriousBug("Failed to create GribHandle from Field", Here());

    return clone(fh,*gh);
}

void Grib::clone( const FieldSet& fields, const PathName& src, const PathName& opath )
{
    bool overwrite = true;

    if( opath.exists() )
        opath.unlink();

    eckit::ScopedPtr<DataHandle> of( opath.fileHandle(overwrite) ); AutoClose of_close(*of);

    ASSERT(of);

    of->openForWrite(0);

    for( size_t i = 0; i < fields.size(); ++i )
    {
		Grib::clone(fields[i], src, *of);
    }
}

void Grib::write(const Field& f, const PathName& opath)
{
    FileHandle fh( opath );

    Length len;
    fh.openForWrite(len);

    write(f, fh);

    fh.close();
}

void Grib::clone(const Field& field, const PathName& gridsec, DataHandle& out )
{
    FILE* fh = ::fopen( gridsec.asString().c_str(), "r" );
    if( fh == 0 )
        throw ReadError( std::string("error opening file ") + gridsec );

    int err = 0;
    grib_handle* clone_h = grib_handle_new_from_file(0,fh,&err);
    if( clone_h == 0 || err != 0 )
        throw ReadError( std::string("error reading grib file ") + gridsec );

    GribHandle ch(clone_h);

	GribHandle::Ptr h = Grib::clone( field, ch );

    //    GRIB_CHECK( grib_write_message(h->raw(),fname.asString().c_str(),"w"), 0 );

    // dump the handle to the DataHandle
    const void* buffer = NULL;
    size_t size = 0;

    GRIB_CHECK( grib_get_message( h->raw(), &buffer, &size), 0);

    out.write(buffer, size);
}

GribHandle::Ptr Grib::clone(const Field& f, GribHandle& gridsec )
{
    const size_t npts = f.size();

    int err=0;
    int what = GRIB_SECTION_GRID;

	GribHandle& meta = f.grib();

    grib_handle* h = grib_util_sections_copy( gridsec.raw(), meta.raw(), what, &err);
    GRIB_CHECK(err,"grib_util_sections_copy()");
    ASSERT( h );

    GribHandle::Ptr gh ( new GribHandle( h ) );

    ASSERT( gh );

    gh->setDataValues(f.data<double>(),npts);

	return gh;
}

struct gridspec_to_grib
{
	gridspec_to_grib( const GridSpec& gspec, GribHandle& gh) :
		gspec_(gspec),
		gh_(gh)
	{}

	GribHandle& gh_;
	const GridSpec& gspec_;

	template <typename T>
	void set( std::string spec, std::string grib )
	{
		if( gspec_.has(spec) )
		{
			GribMutator<T>(grib).set( gh_, gspec_[spec] );
		}
	}
};

void Grib::write_gridspec_to_grib(const GridSpec& gspec, GribHandle& gh)
{
	gridspec_to_grib gspec2grib(gspec,gh);

	if (gh.edition() == 1) {
	   // Clear vertical co-ordinates levels, which for GRIB1 are in geometry section
	   grib_set_long(gh.raw(),"PVPresent", 0 );
	   grib_set_long(gh.raw(),"NV", 0 );
	}

	gspec2grib.set<long>( "Ni", "Ni" );
	gspec2grib.set<long>( "Nj", "Nj" );

	gspec2grib.set<long>( "GaussN", "numberOfParallelsBetweenAPoleAndTheEquator" );

	gspec2grib.set<double>( "grib_bbox_n", "latitudeOfFirstGridPointInDegrees" );
	gspec2grib.set<double>( "grid_bbox_s", "latitudeOfLastGridPointInDegrees" );
	gspec2grib.set<double>( "grid_bbox_w", "longitudeOfFirstGridPointInDegrees" );
	gspec2grib.set<double>( "grid_bbox_e", "longitudeOfLastGridPointInDegrees" );

	gspec2grib.set<double>( "grid_lat_inc", "jDirectionIncrementInDegrees" );
	gspec2grib.set<double>( "grid_lon_inc", "iDirectionIncrementInDegrees" );


	gspec2grib.set<double>( "SouthPoleLat", "latitudeOfSouthernPoleInDegrees" );
	gspec2grib.set<double>( "SouthPoleLon", "longitudeOfSouthernPoleInDegrees" );
	gspec2grib.set<double>( "SouthPoleRotAngle", "angleOfRotation" );

	// Polar stereo Graphic
   gspec2grib.set<long>( "resolutionAndComponentFlag", "resolutionAndComponentFlag" );
   gspec2grib.set<long>( "Nx", "Nx" );
   gspec2grib.set<long>( "Ny", "Ny" );
   gspec2grib.set<long>( "numberOfDataPoints", "numberOfDataPoints" );

   gspec2grib.set<double>( "Dx", "DxInMetres" );
   gspec2grib.set<double>( "Dy", "DyInMetres" );
   gspec2grib.set<double>( "orientationOfTheGrid", "orientationOfTheGridInDegrees" );
   gspec2grib.set<double>( "LaD", "LaDInDegrees" );
   gspec2grib.set<double>( "latitudeOfFirstGridPoint", "latitudeOfFirstGridPointInDegrees" );
   gspec2grib.set<double>( "longitudeOfFirstGridPoint", "longitudeOfFirstGridPointInDegrees" );

   gspec2grib.set<bool>( "iScansPositively", "iScansPositively" );
   gspec2grib.set<bool>( "jScansPositively", "jScansPositively" );
   gspec2grib.set<bool>( "southPoleOnProjectionPlane", "southPoleOnProjectionPlane" );
}

//------------------------------------------------------------------------------------------------------

} // namespace io
} // namespace atlas

#else
#warning "Missing eckit::grib so cannot generate grib output"
#endif
