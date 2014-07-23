/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
#include "atlas/grid/StackGribFile.h"
#include "eckit/log/Log.h"

using namespace std;
using namespace eckit;

//#define DEBUG 1

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

StackGribFile::StackGribFile(const std::string& path) : theGribFile_(PathName(path)), handle_(0)
{
   init(path);
}

StackGribFile::StackGribFile(const eckit::PathName& pathname) : theGribFile_(pathname), handle_(0)
{
   init(pathname.asString());
}

void StackGribFile::init(const std::string& path)
{
#ifdef DEBUG
   Log::info() << "StackGribFile::init Create a grib handle" << std::endl;
#endif
   int err;
   handle_ = grib_handle_new_from_file(0,theGribFile_,&err);
   if (err != 0) Log::error() <<  "StackGribFile::init: grib_handle_new_from_file error " << err << " for file " << path << endl;
}

StackGribFile::~StackGribFile()
{
#ifdef DEBUG
   Log::info()  << "StackGribFile::~StackGribFile() close the grib file and delete handle" << std::endl;
#endif
   int err = grib_handle_delete(handle_);
   if (err != 0) Log::error() <<  "StackGribFile::~StackGribFile(): grib_handle_delete failed " << err << endl;
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
