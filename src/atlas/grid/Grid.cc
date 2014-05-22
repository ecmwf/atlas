/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/Mesh.hpp"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Tesselation.h"

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace atlas
