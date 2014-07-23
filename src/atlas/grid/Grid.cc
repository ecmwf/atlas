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

#include "atlas/mesh/Mesh.hpp"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Tesselation.h"

using namespace eckit;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

Grid::Ptr Grid::create(const Params& p )
{
	return Grid::Ptr( Factory<Grid>::instance().get( p["gridType"] ).create(p) );
}

Grid::Ptr Grid::create(const GridSpec&)
{
	NOTIMP;
}

Grid::Grid()
{
}

Grid::~Grid()
{
}

Mesh& Grid::mesh()
{
    if( !mesh_ )
    {
        mesh_.reset( new Mesh );
        Tesselation::build_mesh( *this, *mesh_ );
    }

    return *mesh_;
}

const Mesh& Grid::mesh() const
{
     if( !mesh_ )
     {
         mesh_.reset( new Mesh );
         Tesselation::build_mesh( *this, *mesh_ );
     }

     return *mesh_;
}

//------------------------------------------------------------------------------------------------------

Grid::BoundBox Grid::makeGlobalBBox()
{
	return Grid::BoundBox( 90., -90., 359.9999999, 0. );
}

Grid::BoundBox Grid::makeBBox(const Params& p)
{
	return BoundBox( p["area_n"],
					 p["area_s"],
					 p["area_e"],
					 p["area_w"] );
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas
