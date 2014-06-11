/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_GribWrite_h
#define atlas_GribWrite_h

#include "atlas/mesh/Mesh.hpp"
#include "atlas/grid/FieldSet.h"

//------------------------------------------------------------------------------------------------------

namespace eckit { class GribHandle; }
namespace eckit { class PathName; }

namespace atlas {

//------------------------------------------------------------------------------------------------------

class GribWrite {

public: // methods

    static void write( const atlas::grid::FieldSet& field, const eckit::PathName& opath  );

    static void clone( const atlas::grid::FieldSet& field, const eckit::PathName& src, const eckit::PathName& opath  );

private: // methods

    static void clone( const atlas::grid::FieldHandle& field, const eckit::PathName& src, const eckit::PathName& opath );

    static eckit::GribHandle* clone(const grid::FieldHandle &field, eckit::GribHandle& source );

};

//---------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif

