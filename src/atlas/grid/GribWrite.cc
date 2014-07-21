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
#include "eckit/config/Resource.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/grib/GribHandle.h"
#include "eckit/utils/Translator.h"
#include "eckit/memory/ScopedPtr.h"
#include "eckit/io/DataHandle.h"
#include "eckit/io/FileHandle.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/parser/StringTools.h"

#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/Parameters.hpp"

#include "atlas/grid/FieldSet.h"
#include "atlas/grid/GribWrite.h"
#include "atlas/grid/GridSpec.h"
#include "atlas/grid/StackGribFile.h"

//------------------------------------------------------------------------------------------------------

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::grid;

namespace atlas {

static std::string map_short_name_to_grib_sample_file(const std::string& short_name, long editionNumber);

//------------------------------------------------------------------------------------------------------

GribHandle::Ptr GribWrite::create_handle( const Grid& grid, long edition )
{
    // determine choice of editionNumber from a resorce

    if( edition == 0 )
        edition = Resource<unsigned>( "NewGribEditionNumber", 2 );

    // From the Grid get the Grid Spec
    eckit::ScopedPtr< GridSpec > gridspec( grid.spec() );

    grib_handle* gh = 0;
    std::string sample_file;

    // first match GridSpec short names, directly to a samples file

    sample_file = map_short_name_to_grib_sample_file( gridspec->short_name(), edition );
    if( !sample_file.empty() )
    {
        gh = grib_handle_new_from_samples(0,sample_file.c_str() );
        if( !gh )
            throw SeriousBug( "Failed to create GribHandle from sample: " + sample_file, Here() );
    }
    else // if this fails, then try looking on disk
    {
        sample_file = GribWrite::grib_sample_file(*gridspec,edition);
        if (!sample_file.empty())
        {
            gh = grib_handle_new_from_samples(0,sample_file.c_str() );
            if( !gh )
                throw SeriousBug( "Failed to create GribHandle from sample: " + sample_file, Here() );
        }
    }

    if(!gh)
        throw SeriousBug( "Failed to create GribHandle from Grib", Here() );

    return GribHandle::Ptr( new GribHandle(gh) );
}


void GribWrite::determine_grib_samples_dir(std::vector<std::string>& sample_paths)
{
   char* the_paths = grib_samples_path(NULL);
   if (the_paths) {
      // Expect <path1>:<path2>:<path3:
      // TODO: Need abstraction for path separator.
      sample_paths = StringTools::split(":", std::string(the_paths));
      return;
   }

   char* include_dir = getenv("GRIB_API_INCLUDE");
   if (!include_dir) throw SeriousBug(string("grib_samples_path(NULL) returned a NULL path"),Here()) ;

   std::string grib_include_dir(include_dir);
   if (grib_include_dir.find("grib_api") == std::string::npos) {
      // "grib-api not found on directory " << grib_include_dir
      return throw SeriousBug(string("grib_samples_path(NULL) returned a NULL path"),Here()) ;
   }

   if (grib_include_dir.find("-I") != std::string::npos) {
      //std::cout << "GRIB_API_INCLUDE=" << grib_include_dir << "\n";
      grib_include_dir.erase(grib_include_dir.begin(),grib_include_dir.begin()+2); // remove -I
   }

   // Handle multiple include dirs
   // If there are any spaces in the string, only take the first include
   size_t space_pos = grib_include_dir.find(" ");
   if (space_pos != std::string::npos) {
      grib_include_dir = grib_include_dir.substr(0,space_pos);
      //std::cout << "GRIB_API_INCLUDE=" << grib_include_dir << "\n";
   }

   // Remove the 'include' and replace with, 'share/grib_api/samples'
   size_t pos = grib_include_dir.find("/include");
   if ( pos == string::npos) {
      // include not found in directory " << grib_include_dir);
      throw SeriousBug(string("grib_samples_path(NULL) returned a NULL path"),Here()) ;
   }

   std::string grib_samples_dir = grib_include_dir.replace(pos,grib_include_dir.length(),"/share/grib_api/samples");
   //std::cout << " GRIB SAMPLES=" << grib_include_dir << "\n";
   sample_paths.push_back( grib_samples_dir );
}

bool match_grid_spec_with_sample_file(
         const GridSpec& the_grid_spec,
         grib_handle* handle,
         long editionNumber,
         const std::string& file_path)
{
   char string_value[64];
   size_t len = sizeof(string_value)/sizeof(char);
   int err = grib_get_string(handle,"gridType",string_value,&len);
   if (err != 0) {
      //Log::error() << "GribWrite::match_grid_spec_with_sample_file, grib_get_string(gridType) failed for \nfile " << file_path << " IGNORING !! " << std::endl;
      return false;
   }
   std::string grib_grid_type = string_value;
   if ( the_grid_spec.grid_type() != grib_grid_type ) {
      //Log::info() << "grid_type in GridSpec " << the_grid_spec.grid_type() << " does not match " << grib_grid_type << " in samples file " << file_path << " IGNORING " << std::endl;
      return false;
   }

   eckit::Value spec_nj = the_grid_spec.get("Nj");
   if (!spec_nj.isNil()) {
      long the_spec_nj = spec_nj;
      long grib_nj = 0;
      if (grib_get_long(handle,"Nj",&grib_nj) == 0 ) {
         if (the_spec_nj != grib_nj ) {
            //Log::info() << "GribWrite::match_grid_spec_with_sample_file, Nj in GridSpec " << the_spec_nj << " does not match  " << grib_nj << " in samples file " << file_path << " IGNORING " << std::endl;
            return false;
         }
      }
   }
   eckit::Value spec_ni = the_grid_spec.get("Ni");
   if (!spec_ni.isNil()) {
      long the_spec_ni = spec_ni;
      long grib_ni = 0;
      if (grib_get_long(handle,"Ni",&grib_ni) == 0 ) {
         if (the_spec_ni != grib_ni ) {
            //Log::info() << "GribWrite::match_grid_spec_with_sample_file, Ni in GridSpec " << the_spec_ni << " does not match  " << grib_ni << " in samples file " << file_path << " IGNORING " << std::endl;
            return false;
         }
      }
   }

   long grib_editionNumber = 0;
   GRIB_CHECK(grib_get_long(handle,"editionNumber",&grib_editionNumber),0);
   if (grib_editionNumber != editionNumber ) {
      //Log::info() << "GribWrite::match_grid_spec_with_sample_file, the_edition_number passed in " << editionNumber << " does not match grib" << editionNumber << " in samples file " << file_path << " IGNORING " << std::endl;
      return false;
   }

   return true;
}

std::string GribWrite::grib_sample_file( const grid::GridSpec& the_grid_spec, long editionNumber )
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

   // From the grid spec, we will look at the grid samples, and find the closest match
   std::vector<std::string> sample_paths;
   determine_grib_samples_dir(sample_paths);
   if ( sample_paths.empty() ) {
      throw SeriousBug(string("Error no sample paths found"),Here()) ;
   }

   for(size_t path = 0; path < sample_paths.size(); ++path) {

      std::string grib_samples_dir = sample_paths[path];
      if (grib_samples_dir.empty()) {
         throw SeriousBug(string("Error, empty samples path. Could not create handle from grid"),Here()) ;
      }
      PathName dir_path(grib_samples_dir);
      if (!dir_path.exists()) continue;
      if (!dir_path.isDir())  continue;

      std::vector<PathName> files;
      std::vector<PathName> directories;
      dir_path.children(files,directories);
      for(size_t i = 0; i < files.size(); i++) {
         try {
            StackGribFile the_grib_file(std::string(files[i].localPath()));

            std::string grib_sample_file_tmpl = files[i].localPath();
            if (match_grid_spec_with_sample_file(the_grid_spec,the_grib_file.handle(),editionNumber,grib_sample_file_tmpl)) {
               // remove .tmpl extension
               eckit::LocalPathName path(grib_sample_file_tmpl);
               LocalPathName the_base_name = path.baseName(false);
               std::string grib_sample_file = the_base_name.localPath();
               return grib_sample_file;
            }
         }
         catch ( const std::exception & ex ) {
            Log::info() << files[i].localPath() << " " << ex.what() << std::endl;
         }
      }
   }

   Log::info() << "Could find grib samples match for grid_spec " << the_grid_spec << std::endl;
   return std::string();
}

static std::string map_short_name_to_grib_sample_file(const std::string& short_name,long editionNumber)
{
   std::stringstream ss; ss << short_name << "_" << editionNumber;
   std::string the_short_name = ss.str();

   // Short cut, for mapping short name to grib samples file.
   std::map<std::string,std::string> short_name_to_samples_map;
   short_name_to_samples_map.insert( std::make_pair(std::string("QG32_1"),std::string("reduced_gg_pl_32_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG32_2"),std::string("reduced_gg_pl_32_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG48_1"),std::string("reduced_gg_pl_48_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG48_2"),std::string("reduced_gg_pl_48_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG80_1"),std::string("reduced_gg_pl_80_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG80_2"),std::string("reduced_gg_pl_80_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG128_1"),std::string("reduced_gg_pl_128_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG128_2"),std::string("reduced_gg_pl_128_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG160_1"),std::string("reduced_gg_pl_160_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG160_2"),std::string("reduced_gg_pl_160_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG200_1"),std::string("reduced_gg_pl_200_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG200_2"),std::string("reduced_gg_pl_200_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG256_1"),std::string("reduced_gg_pl_256_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG256_2"),std::string("reduced_gg_pl_256_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG320_1"),std::string("reduced_gg_pl_320_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG320_2"),std::string("reduced_gg_pl_320_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG400_1"),std::string("reduced_gg_pl_400_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG400_2"),std::string("reduced_gg_pl_400_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG512_1"),std::string("reduced_gg_pl_512_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG512_2"),std::string("reduced_gg_pl_512_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG640_1"),std::string("reduced_gg_pl_640_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG640_2"),std::string("reduced_gg_pl_640_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG1024_1"),std::string("reduced_gg_pl_1024_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG1024_2"),std::string("reduced_gg_pl_1024_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG1280_1"),std::string("reduced_gg_pl_1280_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG1280_2"),std::string("reduced_gg_pl_1280_grib2")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG2000_1"),std::string("reduced_gg_pl_2000_grib1")) );
   short_name_to_samples_map.insert( std::make_pair(std::string("QG2000_2"),std::string("reduced_gg_pl_2000_grib2")) );

   std::map<std::string,std::string>::const_iterator i = short_name_to_samples_map.find(the_short_name);
   if (i != short_name_to_samples_map.end()) {
      return (*i).second;
   }
   return std::string();
}



void GribWrite::write( const FieldSet& fields, const PathName& opath )
{
    for( size_t i = 0; i < fields.size(); ++i )
    {
        PathName pi( opath.asString() + "." + Translator<size_t,std::string>()(i) );
        GribWrite::write(fields[i], pi);
    }
}

void GribWrite::write(const FieldHandle& fh, DataHandle& dh)
{
    GribHandle::Ptr gh = GribWrite::create_handle( fh.grid(), fh.grib().edition() );

    if( !gh )
        throw SeriousBug("Failed to create GribHandle from FieldHandle", Here());

    if( !gh->raw() )
        throw SeriousBug("Failed to create GribHandle from FieldHandle", Here());

    GribHandle::Ptr h = clone(fh,*gh);

    // dump the handle to the DataHandle
    const void* buffer = NULL;
    size_t size = 0;

    GRIB_CHECK( grib_get_message( h->raw(), &buffer, &size), 0);

    dh.write(buffer, size);
}

GribHandle::Ptr GribWrite::write(const FieldHandle& fh)
{
    GribHandle::Ptr gh = GribWrite::create_handle( fh.grid(), fh.grib().edition() );

    if( !gh )
        throw SeriousBug("Failed to create GribHandle from FieldHandle", Here());

    return clone(fh,*gh);
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

void GribWrite::write(const FieldHandle& f, const PathName& opath)
{
    FileHandle fh( opath );

    Length len;
    fh.openForWrite(len);

    write(f, fh);

    fh.close();
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

    GribHandle::Ptr h = GribWrite::clone( field, ch );

    //    GRIB_CHECK( grib_write_message(h->raw(),fname.asString().c_str(),"w"), 0 );

    // dump the handle to the DataHandle
    const void* buffer = NULL;
    size_t size = 0;

    GRIB_CHECK( grib_get_message( h->raw(), &buffer, &size), 0);

    out.write(buffer, size);
}

GribHandle::Ptr GribWrite::clone(const FieldHandle& field, GribHandle& gridsec )
{
    const Field& f = field.data();
    const size_t npts = f.size();

    // check number of points matches

    size_t nb_nodes = gridsec.getNbDataPoints();
    ASSERT( nb_nodes == f.size() );

    ///@todo move this to the eckit::grib interface
    int err=0;
    int what = GRIB_SECTION_GRID;

    GribHandle& meta = field.grib();

    grib_handle* h = grib_util_sections_copy( gridsec.raw(), meta.raw(), what, &err);
    GRIB_CHECK(err,"grib_util_sections_copy()");
    ASSERT( h );

    GribHandle::Ptr gh ( new GribHandle( h ) );

    ASSERT( gh );

    gh->setDataValues(f.data<double>(),npts);

    return gh;
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

