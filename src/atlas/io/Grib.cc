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

/// @todo Remove these once we use only eckit::grib
// Do not add following file because not present in grib_api 1.13
// --> #include "grib_api_config.h"
#include "grib_api.h"

#include "eckit/eckit_config.h"
#include "eckit/config/Resource.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/utils/Translator.h"
#include "eckit/memory/ScopedPtr.h"
#include "eckit/io/DataHandle.h"
#include "eckit/io/FileHandle.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/parser/StringTools.h"

#include "eckit/grib/GribParams.h"
#include "eckit/grib/GribHandle.h"
#include "eckit/grib/GribMutator.h"
#include "eckit/grib/GribAccessor.h"

#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Parameters.h"

#include "atlas/FieldSet.h"

#include "atlas/io/Grib.h"
#include "atlas/io/PointCloud.h"

#include "atlas/grids/Unstructured.h"
#include "atlas/grids/ReducedGaussianGrid.h"
#include "atlas/grids/GaussianGrid.h"
#include "atlas/grids/LonLatGrid.h"
#include "atlas/grids/ReducedLonLatGrid.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace eckit::grib;
using namespace atlas;

namespace atlas {
namespace io {

static std::string match_grid_spec_with_sample_file( const eckit::Properties& g_spec, long edition );

//------------------------------------------------------------------------------------------------------

Grid::Ptr make_grid( const PathName& path )
{
  // read start of file (using buffer) to determine file format
  Buffer buff(64);
  DataHandle* dh = path.fileHandle();
  dh->openForRead();
  size_t len = dh->read(buff,buff.size());
  ASSERT(len);
  dh->close();

  if ((len>=4) && (0==strncmp(buff,"GRIB",4)))
  {
    FILE* fh = ::fopen( path.asString().c_str(), "r" );
    if( fh == 0 )
        throw ReadError( std::string("error opening file ") + path.asString() );

    int err = 0;
    grib_handle* h;

    h = grib_handle_new_from_file(0,fh,&err);

    if( h == 0 || err != 0 )
        throw ReadError( std::string("error reading grib file ") + path.asString() );

    if( ::fclose(fh) == -1 )
        throw ReadError( std::string("error closing file ") + path.asString() );

    GribHandle gh(h);
    Grid::Ptr g ( Grib::create_grid( gh ) );
    ASSERT( g );

    return g;
  }

  // attempt to read PointCloud format

  if ((len>=10) && (0==strncmp(buff,"PointCloud",10)))
  {
    Mesh::Ptr mesh( io::PointCloud::read(path) );
    Grid::Ptr g( new grids::Unstructured(*mesh) );
    return g;
  }

  // Missing return... fail here
  Grid::Ptr g;
  ASSERT( g );
  return g;
}

Grid::Ptr Grib::create_grid(GribHandle& gh)
{
	eckit::Properties* gp = GribParams::create(gh);

	ASSERT( gp );

  Config params(*gp);
	Grid::Ptr grid( Grid::create( params ) );

  delete gp;
  return grid;
}

GribHandle::Ptr Grib::create_handle( const Grid& grid, long edition )
{
  // determine choice of editionNumber from a resource
  if( edition == 0 )
    edition = Resource<unsigned>( "NewGribEditionNumber", 2 );


  // From the Grid get the Grid Spec
  eckit::Properties grid_spec = grid.spec();


  grib_handle* gh = 0;

  // Try to match grid with samples file *WITHOUT* going to disk
  // This should succeed most of the time.
  // Will fail for Polar Stereographic, since the sample for these are not correct
  std::string sample_file = match_grid_spec_with_sample_file(grid_spec,edition);

  if( !gh && !sample_file.empty() )
  {
    gh = grib_handle_new_from_samples(0,sample_file.c_str() );

    if( !gh )
      throw SeriousBug( "Failed to create GribHandle from sample: " + sample_file, Here() );
  }

  /// scan SAMPLE files directories and match to a grib

  if( !gh )
  {
    sample_file = Grib::grib_sample_file(grid_spec,edition);

    if( sample_file.empty() )
      throw BadParameter("Could not find GRIB sample for grid " + grid.shortName(), Here() );

    gh = grib_handle_new_from_samples(0,sample_file.c_str() );

    if( !gh )
      throw SeriousBug( "Failed to create GribHandle from sample: " + sample_file, Here() );
  }

  // finally create handle

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

#if ( GRIB_API_VERSION >= 11300 )
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
         path += "/share/grib_api/samples";
         sample_paths.push_back(path);
         return;
      }
	   throw SeriousBug("GRIB_SAMPLES_PATH and GRIB_API_PATH not defined, please call: module load grib_api",Here()) ;
   }

   sample_paths = StringTools::split(":", std::string(samples_dir));
}


static std::string match_grid_spec_with_sample_file( const eckit::Properties& g_spec, long edition )
{
   // First get the grid_type
   std::string grid_type = g_spec["grid_type"];

   // For regular gaussian grids
   if (grid_type == grids::GaussianGrid::grid_type_str() ) {


      // regular_gg_ml_grib1.tmpl  --> eckit::Properties[ (GaussN,32)(Ni,128)(Nj,64)(uid,regular_gg_32) ]
      // regular_gg_ml_grib2.tmpl  --> eckit::Properties[ (GaussN,32)(Ni,128)(Nj,64)(uid,regular_gg_32) ]
      // regular_gg_pl_grib1.tmpl  --> eckit::Properties[ (GaussN,32)(Ni,128)(Nj,64)(uid,regular_gg_32) ]
      // regular_gg_pl_grib2.tmpl  --> eckit::Properties[ (GaussN,32)(Ni,128)(Nj,64)(uid,regular_gg_32) ]
      // regular_gg_sfc_grib1.tmpl --> eckit::Properties[ (GaussN,32)(Ni,128)(Nj,64)(uid,regular_gg_32) ]
      // regular_gg_sfc_grib2.tmpl --> eckit::Properties[ (GaussN,32)(Ni,128)(Nj,64)(uid,regular_gg_32) ]
      if (edition == 1) return "regular_gg_ml_grib1";
      return "regular_gg_ml_grib2";
   }


   // For reduced lat long, Choice of two only, if we dont have 501 points number of pts north to south, were out of luck.
   if (grid_type == grids::ReducedLonLatGrid::grid_type_str() ) {
      if (edition == 1) return "reduced_ll_sfc_grib1";
      return "reduced_ll_sfc_grib2";
   }

   if (grid_type == grids::LonLatGrid::grid_type_str() ) {

      // regular_ll_pl_grib1.tmpl  --> eckit::Properties[ (Ni,16)(Nj,31) (lat_inc,2)(lon_inc,2)(uid,regular_ll_31_16) ]
      // regular_ll_pl_grib2.tmpl  --> eckit::Properties[ (Ni,16)(Nj,31) (lat_inc,2)(lon_inc,2)(uid,regular_ll_31_16) ]
      // regular_ll_sfc_grib1.tmpl --> eckit::Properties[ (Ni,16)(Nj,31) (lat_inc,2)(lon_inc,2)(uid,regular_ll_31_16) ]
      // regular_ll_sfc_grib2.tmpl --> eckit::Properties[ (Ni,16)(Nj,31) (lat_inc,2)(lon_inc,2)(uid,regular_ll_31_16) ]
      if (edition == 1) return "regular_ll_pl_grib1";
      return "regular_ll_pl_grib2";;
   }

   /// NOTE: rotated_ll doesn't exist anymore in Atlas as a standalone grid

   if (grid_type == "rotated_ll" ) {

      // rotated_ll_pl_grib1.tmpl  --> eckit::Properties[ (Ni,16)(Nj,31)(SouthPoleLat,0)(SouthPoleLon,0)(SouthPoleRotAngle,0)(lat_inc,2)(lon_inc,2)(uid,rotated_ll_31) ]
      // rotated_ll_pl_grib2.tmpl  --> eckit::Properties[ (Ni,16)(Nj,31)(SouthPoleLat,0)(SouthPoleLon,0)(SouthPoleRotAngle,0)(lat_inc,2)(lon_inc,2)(uid,rotated_ll_31) ]
      // rotated_ll_sfc_grib1.tmpl --> eckit::Properties[ (Ni,16)(Nj,31)(SouthPoleLat,0)(SouthPoleLon,0)(SouthPoleRotAngle,0)(lat_inc,2)(lon_inc,2)(uid,rotated_ll_31) ]
      // rotated_ll_sfc_grib2.tmpl --> eckit::Properties[ (Ni,16)(Nj,31)(SouthPoleLat,0)(SouthPoleLon,0)(SouthPoleRotAngle,0)(lat_inc,2)(lon_inc,2)(uid,rotated_ll_31) ]
      if (edition == 1) return "rotated_ll_pl_grib1";
      return "rotated_ll_pl_grib2";
   }

   // Cases we might fail to match:
   //  * GRIB samples for polar stereographic are in-correct,
   //  * Can't handle spherical harmonics yet

   return std::string(); // returning emtpy string when match fails
}

bool check_grid_spec_matches_sample_file( const eckit::Properties& g_spec, long edition,	const eckit::PathName& fpath )
{
	GribHandle gh( fpath );

   if ( g_spec["grid_type"] != gh.gridType() )
   {
      //Log::info() << "grid_type in eckit::Properties " << g_spec.gridType() << " does not match " << grib_grid_type << " in samples file " << file_path << " IGNORING " << std::endl;
      return false;
   }

   if( gh.edition() != edition )
   {
      //Log::info() << "Grib::match_grid_spec_with_sample_file, edition_number passed in " << edition << " does not match grib" << edition << " in samples file " << file_path << " IGNORING " << std::endl;
      return false;
   }

   if (g_spec["grid_type"] == grids::ReducedGaussianGrid::grid_type_str() )
   {
      if (g_spec.has("N"))
      {
         long grid_gausn = g_spec.get("N");
         long grib_gausn = GribAccessor<long>("numberOfParallelsBetweenAPoleAndTheEquator")(gh);

         if (grid_gausn != grib_gausn)
            return false;

      }
   }

   return true;
}

std::string Grib::grib_sample_file( const eckit::Properties& g_spec, long edition )
{
    // Note: many of the grib samples files are not UNIQUE in their grid specification:
    // i.e
    //   GRIB2.tmpl                        -> eckit::Properties[ regular_ll, LL31_16_2, Ni:16, Nj:31, typeOfLevel:surface ]
    //   regular_ll_pl_grib2.tmpl          -> eckit::Properties[ regular_ll, LL31_16_2, Ni:16, Nj:31 ]
    //   regular_ll_sfc_grib2.tmpl         -> eckit::Properties[ regular_ll, LL31_16_2, Ni:16, Nj:31 ]
    //
    //   reduced_gg_ml_grib1               -> eckit::Properties[ reduced_gg, QG32_1, Nj:64 ]
    //   reduced_gg_pl_32_grib1            -> eckit::Properties[ reduced_gg, QG32_1, Nj:64 ]
    //   reduced_gg_ml_grib2               -> eckit::Properties[ reduced_gg, QG32_2, Nj:64 ]
    //   reduced_gg_pl_32_grib2            -> eckit::Properties[ reduced_gg, QG32_2, Nj:64 ]
    //
    // Others are just plain wrong, i.e
    //   polar_stereographic_pl_grib2.tmpl -> eckit::Properties[ rotated_ll, RL31_2, Ni:16, Nj:31, editionNumber:2 ]
    //
    // From the grid spec, we will look at the grid samples, and find the closest match

    std::vector<std::string> sample_paths;
    determine_grib_samples_dir(sample_paths);

    if ( sample_paths.empty() )
        throw SeriousBug("Error no sample paths found",Here()) ;

    for(size_t path = 0; path < sample_paths.size(); ++path) {
        std::string grib_samples_dir = sample_paths[path];

        DEBUG_VAR( grib_samples_dir );

        if (grib_samples_dir.empty())
            throw SeriousBug("Error, empty samples path. Could not create handle from grid",Here()) ;

        PathName dir_path(grib_samples_dir);

        if( !dir_path.exists() ) continue;
        if( !dir_path.isDir()  ) continue;

        std::vector<PathName> files;
        std::vector<PathName> directories;
        dir_path.children(files,directories);

        for(size_t i = 0; i < files.size(); i++) {
            // FIXME: we only consider files with .tmpl extension so we don't pick
            // up non-GRIB files in the build directory (Makefile.am etc.)
            const std::string& file = files[i].path();
            if(file.substr(file.find_last_of('.') + 1) == "tmpl" &&
                    check_grid_spec_matches_sample_file(g_spec, edition, files[i]))
                return files[i].baseName(false); // remove .tmpl extension
        }
    }

//    Log::warning() << "Could *not* find grib samples match for grid_spec " << g_spec << std::endl;

    return std::string();
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
  GribHandle::Ptr gh;
   NOTIMP; // Can no longer access grid from field directly
   // ---> gh = Grib::create_handle( fh.grid() );

   if( !gh )
      throw SeriousBug("Failed to create GribHandle from Field", Here());

   if( !gh->raw() )
      throw SeriousBug("Failed to create GribHandle from Field", Here());

   GribHandle::Ptr h ( clone(fh,*gh) );

   h->write(dh);
 }

GribHandle::Ptr Grib::write(const Field& fh)
{
  NOTIMP;
  return GribHandle::Ptr();


//  long edition = fh.grib() ? fh.grib()->edition() : 2;

//	GribHandle::Ptr gh = Grib::create_handle( fh.grid(), edition );

//	if( !gh )
//		throw SeriousBug("Failed to create GribHandle from Field", Here());

//	return GribHandle::Ptr( clone(fh,*gh) );
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

GribHandle* Grib::copy_metadata(const GribHandle& from, GribHandle& to)
{
	int err=0;
	int what = GRIB_SECTION_PRODUCT | GRIB_SECTION_LOCAL;

	grib_handle* h = grib_util_sections_copy( from.raw(), to.raw(), what, &err);
	GRIB_CHECK(err,"grib_util_sections_copy()");
	ASSERT( h );

	return new GribHandle( h );
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

   GribHandle::Ptr h ( Grib::clone( field, ch ) );

   h->write(out);
}

GribHandle* Grib::clone(const Field& f, GribHandle& gridsec )
{
    GribHandle* gh = 0;

    NOTIMP;
    return gh;

//    if( f.grib() )
//    {
//      const GribHandle& meta = *(f.grib());
//      int err=0;
//      int what = GRIB_SECTION_GRID;

//      grib_handle* h = grib_util_sections_copy( gridsec.raw(), meta.raw(), what, &err);
//      GRIB_CHECK(err,"grib_util_sections_copy()");
//      ASSERT( h );
//      gh = new GribHandle( h );
//    }
//    else
//    {
//      gh = gridsec.clone();
//    }

//    ASSERT( gh );

//    gh->setDataValues(f.data<double>(),npts);

//    return gh;
}

struct gridspec_to_grib
{
  gridspec_to_grib( const eckit::Properties& gspec, GribHandle& gh) :
    gh_(gh),
    gspec_(gspec)
  {}

  GribHandle& gh_;
  const eckit::Properties& gspec_;

  template <typename T>
  void set( std::string spec, std::string grib )
  {
    if( gspec_.has(spec) )
    {
      GribMutator<T>(grib).set( gh_, gspec_[spec] );
    }
  }
};

void Grib::write_gridspec_to_grib(const eckit::Properties& gspec, GribHandle& gh)
{
  gridspec_to_grib gspec2grib(gspec,gh);

  if (gh.edition() == 1) {
    // Clear vertical co-ordinates levels, which for GRIB1 are in geometry section
    grib_set_long(gh.raw(),"PVPresent", 0 );
    grib_set_long(gh.raw(),"NV", 0 );
  }

  gspec2grib.set<long>( "nlon", "Ni" );
  gspec2grib.set<long>( "nlat", "Nj" );

  // N number for (Reduced) Gaussian Grids
  if( gspec["grid_type"] == grids::GaussianGrid::grid_type_str() ||
      gspec["grid_type"] == grids::ReducedGaussianGrid::grid_type_str() )
    gspec2grib.set<long>( "N", "numberOfParallelsBetweenAPoleAndTheEquator" );

  gspec2grib.set<double>( "bbox_n", "latitudeOfFirstGridPointInDegrees" );
  gspec2grib.set<double>( "bbox_s", "latitudeOfLastGridPointInDegrees" );
  gspec2grib.set<double>( "bbox_w", "longitudeOfFirstGridPointInDegrees" );
  gspec2grib.set<double>( "bbox_e", "longitudeOfLastGridPointInDegrees" );

  gspec2grib.set<double>( "lat_inc", "jDirectionIncrementInDegrees" );
  gspec2grib.set<double>( "lon_inc", "iDirectionIncrementInDegrees" );


  gspec2grib.set<double>( "SouthPoleLat", "latitudeOfSouthernPoleInDegrees" );
  gspec2grib.set<double>( "SouthPoleLon", "longitudeOfSouthernPoleInDegrees" );
  gspec2grib.set<double>( "SouthPoleRotAngle", "angleOfRotation" );

  // Polar stereo Graphic
   gspec2grib.set<long>( "resolutionAndComponentFlag", "resolutionAndComponentFlag" );
   gspec2grib.set<long>( "Nx", "Nx" );
   gspec2grib.set<long>( "Ny", "Ny" );
   gspec2grib.set<long>( "npts", "numberOfDataPoints" );

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
 // Missing eckit::grib so GRIB IO is not supported in Atlas
#endif
