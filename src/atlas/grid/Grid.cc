/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/memory/Factory.h"
#include "eckit/memory/Builder.h"
#include "eckit/config/Resource.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Tesselation.h"
#include "atlas/grid/GridSpecParams.h"

using namespace eckit;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

Grid::Ptr Grid::create(const Params& p )
{
	return Grid::Ptr( Factory<Grid>::instance().get( p["grid_type"] ).create(p) );
}


Grid::Ptr Grid::create(const GridSpec& g)
{
	GridSpecParams p( g );
   return Grid::Ptr( Factory<Grid>::instance().get( p["grid_type"] ).create(p) );
}

Grid::Grid()
{
}

Grid::~Grid()
{
}

void Grid::buildMesh() const
{
	if( ! mesh_ )
	{
		mesh_.reset( new Mesh() );
		mesh_->grid( const_cast<Grid&>(*this) );
		Tesselation::build_mesh( *this, *mesh_ );
	}
}

Mesh& Grid::mesh()
{
	buildMesh();
    return *mesh_;
}

const Mesh& Grid::mesh() const
{
	buildMesh();
	return *mesh_;
}

//------------------------------------------------------------------------------------------------------

Grid::BoundBox Grid::makeGlobalBBox()
{
	return Grid::BoundBox( 90., -90., 360. - degrees_eps(), 0. );
}

Grid::BoundBox Grid::makeBBox(const Params& p)
{
	if( p.get("grid_bbox_s").isNil() )
		return Grid::makeGlobalBBox();

	return BoundBox( p["grib_bbox_n"],
					 p["grid_bbox_s"],
					 p["grid_bbox_e"],
			p["grid_bbox_w"] );
}

double Grid::degrees_eps()
{
	/// default is 1E-3 because
	/// some bugs in IFS means we need a lower resolution epsilon when decoding from grib2

	static double eps = eckit::Resource<double>( "$ATLAS_GRID_DEGREES_EPSILON;GridDegreesEps", 1E-3 );
	return eps;
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas
