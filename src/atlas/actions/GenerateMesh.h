/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_GenerateMesh_h
#define atlas_GenerateMesh_h

#include "atlas/mesh/Mesh.h"

namespace atlas {
namespace actions {

Mesh* generate_reduced_gaussian_grid( const std::string& identifier );
Mesh* generate_reduced_gaussian_grid( const std::vector<long>& nlon );
Mesh* generate_full_gaussian_grid( int nlon, int nlat );
Mesh* generate_regular_grid( int nlon, int nlat ); // must be even numbers

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C"
{
  Mesh* atlas__generate_reduced_gaussian_grid (char* identifier);
  Mesh* atlas__generate_full_gaussian_grid (int nlon, int nlat);
  Mesh* atlas__generate_latlon_grid (int nlon, int nlat);
  Mesh* atlas__generate_custom_reduced_gaussian_grid (int nlon[], int nlat);
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

#endif
