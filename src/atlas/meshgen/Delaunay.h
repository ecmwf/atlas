/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_meshgen_Delaunay_h
#define atlas_meshgen_Delaunay_h

#include "atlas/meshgen/MeshGenerator.h"

namespace atlas {

class Mesh;
class Grid;

namespace meshgen {

//------------------------------------------------------------------------------------------------------

class Delaunay : public MeshGenerator {

public:

    Delaunay();

    virtual ~Delaunay();

    virtual void tesselate(const Grid& g, Mesh& mesh) const;

};

//------------------------------------------------------------------------------------------------------

} // namespace meshgen
} // namespace atlas

#endif
