/*
 * (C) Copyright 1996-2012 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "eckit/runtime/Tool.h"
#include "eckit/config/Resource.h"
#include "atlas/atlas.h"

//-----------------------------------------------------------------------------

using namespace eckit;
namespace atlas {

//-----------------------------------------------------------------------------

class Version : public Tool {
public:

    Version(int argc,char **argv): Tool(argc,argv) {}

    ~Version() {}

    virtual void run();
};

std::string print(bool v)
{
  return v ? "ON" : "OFF";
}

void Version::run()
{
  if( Resource<bool>("--version",false) )
  {
    Log::info() << atlas_version() << std::endl;
    return;
  }
  else if( Resource<bool>("--git",false) )
  {
    Log::info() << atlas_git_sha1_abbrev(12) << std::endl;
    return;
  }
  else if( Resource<bool>("--info",false) )
  {
    Log::info() << "atlas version (" << atlas_version_str() << "), "
                << "git-sha1 "<< atlas_git_sha1_abbrev(7) << std::endl;
    Log::info() << std::endl;
    Log::info() << "  Build:" << std::endl;
    Log::info() << "    build type      : " << ATLAS_BUILD_TYPE << std::endl
                << "    timestamp       : " << ATLAS_BUILD_TIMESTAMP << std::endl
                << "    op. system      : " << ATLAS_OS_NAME << " (" << ATLAS_OS_STR << ")"  << std::endl
                << "    processor       : " << ATLAS_SYS_PROCESSOR  << std::endl
                << "    c compiler      : " << ATLAS_C_COMPILER_ID << " " << ATLAS_C_COMPILER_VERSION << std::endl
                << "      flags         : " << ATLAS_C_FLAGS << std::endl
                << "    c++ compiler    : " << ATLAS_CXX_COMPILER_ID << " " << ATLAS_CXX_COMPILER_VERSION << std::endl
                << "      flags         : " << ATLAS_CXX_FLAGS << std::endl
#ifndef EC_HAVE_FORTRAN
                    << "    fortran         : NO " << std::endl
#else
                << "    fortran compiler: " << ATLAS_Fortran_COMPILER_ID << " " << ATLAS_Fortran_COMPILER_VERSION << std::endl
                << "      flags         : " << ATLAS_Fortran_FLAGS << std::endl
#endif
                << std::endl;

    bool feature_fortran(false);
    bool feature_MPI(false);
    bool feature_OpenMP(false);
    bool feature_Grib(false);
    bool feature_Trans(false);
    bool feature_Tesselation(false);
#ifdef ATLAS_HAVE_FORTRAN
      feature_fortran = true;
#endif
#ifdef ATLAS_HAVE_MPI
      feature_MPI = true;
#endif
#ifdef ATLAS_HAVE_OMP
      feature_OpenMP = true;
#endif
#ifdef ECKIT_HAVE_GRIB
      feature_Grib = true;
#endif
#ifdef ATLAS_HAVE_TRANS
      feature_Trans = true;
#endif
#ifdef ATLAS_HAVE_TESSELATION
      feature_Tesselation = true;
#endif
      Log::info() << "  Features:" << std::endl;
    Log::info() << "    Fortran     : " << print(feature_fortran) << std::endl
                << "    MPI         : " << print(feature_MPI) << std::endl
                << "    OpenMP      : " << print(feature_OpenMP) << std::endl
                << "    Grib        : " << print(feature_Grib) << std::endl
                << "    Trans       : " << print(feature_Trans) << std::endl
                << "    Tesselation : " << print(feature_Tesselation) << std::endl
                << "    gidx_t      : " << ATLAS_BITS_GLOBAL << " bit integer" << std::endl
                << std::flush;
    return;
  }
  else if( Resource<bool>("--help",false) )
  {
    Log::info() <<
    "NAME\n"
    "       atlas - Framework for parallel flexible data structures on the sphere\n"
    "\n"
    "SYNOPSIS\n"
    "       atlas [--help] [--version] [--git] [--info]\n"
    "\n"
    "DESCRIPTION\n"
    "       Framework for parallel flexible data structures on the sphere.\n"
    "\n"
    "OPTIONS\n"
    "       --help\n"
    "           Print this help\n"
    "\n"
    "       --version\n"
    "           Print short version string 'MAJOR.MINOR.PATCH' \n"
    "\n"
    "       --info\n"
    "           Print build configuration anad features\n"
    "\n"
    "AUTHOR\n"
    "       Written by Willem Deconinck.\n"
    "\n"
    "ECMWF                        December 2014"
    << std::endl;
    return;
  }
  else
    Log::info() << "usage: atlas [--help] [--version] [--git] [--info]" << std::endl;
}

} // namespace atlas

//-----------------------------------------------------------------------------

using namespace atlas;

int main(int argc,char **argv)
{
    Version app(argc,argv);
    app.start();
    return 0;
}

