/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_actions_BuildXYZField_h
#define atlas_actions_BuildXYZField_h

#include <string>

namespace atlas {

class FunctionSpace;
class Mesh;
class Field;

namespace actions {

Field& build_xyz_field( Mesh&,          const std::string& name = "xyz" );
Field& build_xyz_field( FunctionSpace&, const std::string& name = "xyz" );

} // namespace actions
} // namespace atlas

#endif
