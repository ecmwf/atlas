/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/meshgen/Delaunay.h"

#include "atlas/Tesselation.h"

namespace atlas {
namespace meshgen {

Delaunay::Delaunay()
{
}

Delaunay::~Delaunay() {
}

void Delaunay::tesselate(const Grid &g, Mesh &mesh) const
{
    Tesselation::delaunay_triangulation(mesh);
}

namespace {
static MeshGeneratorBuilder< Delaunay > __delaunay("Delaunay");
}


} // namespace meshgen
} // namespace atlas

