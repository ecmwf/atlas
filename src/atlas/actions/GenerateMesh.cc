/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mpl/MPL.h"
#include "atlas/meshgen/RGG.h"
#include "atlas/meshgen/RGGMeshGenerator.h"
#include "atlas/actions/GenerateMesh.h"
#include "atlas/util/Debug.h"
using namespace atlas::meshgen;

namespace atlas {
namespace actions {

// ------------------------------------------------------------------

Mesh* generate_mesh (const RGG& rgg)
{
  RGGMeshGenerator generate;
  generate.options.set( "nb_parts", MPL::size() );
  generate.options.set( "part"    , MPL::rank() );
  return generate(rgg);
}

// ------------------------------------------------------------------

Mesh* generate_reduced_gaussian_grid( const std::string& identifier )
{
	RGG* rgg = new_reduced_gaussian_grid(identifier);
	Mesh* mesh = generate_mesh(*rgg);
	delete rgg;
	return mesh;
}

// ------------------------------------------------------------------

Mesh* generate_reduced_gaussian_grid( const std::vector<long>& nlon )
{
	RGG* rgg = new_reduced_gaussian_grid(nlon);
	Mesh* mesh = generate_mesh(*rgg);
	delete rgg;
	return mesh;
}

// ------------------------------------------------------------------

Mesh* generate_regular_grid( int nlon, int nlat )
{
	RGG* rgg = new_regular_latlon_grid(nlon,nlat);
	Mesh* mesh = generate_mesh(*rgg);
	delete rgg;
	return mesh;
}

// ------------------------------------------------------------------

Mesh* generate_full_gaussian_grid( int nlon, int nlat )
{
	RGG* rgg = new_regular_gaussian_grid(nlon,nlat);
	Mesh* mesh = generate_mesh(*rgg);
	delete rgg;
	return mesh;
}

// ------------------------------------------------------------------

// C wrapper interfaces to C++ routines
Mesh* atlas__generate_reduced_gaussian_grid (char* identifier)
{
  return generate_reduced_gaussian_grid( std::string(identifier) );
}

Mesh* atlas__generate_full_gaussian_grid ( int nlon, int nlat )
{
  return generate_full_gaussian_grid( nlon, nlat );
}

Mesh* atlas__generate_latlon_grid ( int nlon, int nlat )
{
  return generate_regular_grid( nlon, nlat );
}

Mesh* atlas__generate_custom_reduced_gaussian_grid( int nlon[], int nlat )
{
  std::vector<long> nlon_vector;
  nlon_vector.assign(nlon,nlon+nlat);
  return generate_reduced_gaussian_grid(nlon_vector);
}

Mesh* atlas__generate_mesh (RGG* rgg)
{
  return generate_mesh(*rgg);
}

// ------------------------------------------------------------------


} // namespace actions
} // namespace atlas
