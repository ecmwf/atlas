/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef eckit_GribWrite_h
#define eckit_GribWrite_h

#include "grib_api.h"

#include "atlas/Mesh.hpp"
#include "atlas/grid/Field.h"

//------------------------------------------------------------------------------------------------------

namespace eckit {

//------------------------------------------------------------------------------------------------------

class GribWrite {
public:

    static void write( atlas::grid::FieldH& field, grib_handle* input_h );

    static void clone( atlas::grid::FieldH& field, const std::string& source, const std::string& output );

    static grib_handle* clone( atlas::grid::FieldH& field, grib_handle* source );

};

//---------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif

