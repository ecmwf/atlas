#ifndef atlas_stack_grib_file_H
#define atlas_stack_grib_file_H
/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cstddef>
#include <string>
#include "grib_api.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/io/StdFile.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

// ==========================================================================
// Ensures grib handle is always deleted, in presence of exceptions
class StackGribFile {
public:
   StackGribFile(const std::string& the_file_path);
   StackGribFile(const eckit::PathName& pathname);
   ~StackGribFile();

   grib_handle* handle() const { return handle_;}

private:
   eckit::StdFile theGribFile_;
   grib_handle* handle_;

private:
   void init(const std::string& the_file_path);
   StackGribFile(const StackGribFile&);            // prevent copies
   StackGribFile& operator=(const StackGribFile&); // prevent assignment
};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
