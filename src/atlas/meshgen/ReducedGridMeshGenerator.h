/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef ReducedGridMeshGenerator_h
#define ReducedGridMeshGenerator_h

#include "atlas/Metadata.h"
#include "atlas/meshgen/MeshGenerator.h"

namespace atlas {

class Mesh;
class GridDistribution;

namespace grids { class ReducedGrid; }
namespace meshgen { struct Region; }

}

namespace atlas {
namespace meshgen {

//------------------------------------------------------------------------------------------------------

class ReducedGridMeshGenerator : public MeshGenerator {

public:

  ReducedGridMeshGenerator();

  virtual void tesselate(const Grid& g, Mesh& mesh) const;

  void generate( const grids::ReducedGrid&, const GridDistribution&, Mesh& ) const;
  void generate( const grids::ReducedGrid&, GridDistribution*, Mesh& ) const;
  void generate( const grids::ReducedGrid&, Mesh& ) const;

  Mesh* generate( const grids::ReducedGrid&, const GridDistribution& ) const;
  Mesh* generate( const grids::ReducedGrid&, GridDistribution* ) const;
  Mesh* generate( const grids::ReducedGrid& ) const;

  Mesh* operator()( const grids::ReducedGrid&, const GridDistribution& ) const;
  Mesh* operator()( const grids::ReducedGrid&, GridDistribution* ) const;
  Mesh* operator()( const grids::ReducedGrid& ) const;

private:

  void generate_region( const grids::ReducedGrid&, const std::vector<int>& parts, int mypart, Region& region ) const;

  void generate_mesh( const grids::ReducedGrid&,const std::vector<int>& parts, const Region& region, Mesh& m ) const;

  void generate_global_element_numbering( Mesh& mesh ) const;

public:

  Metadata options;

};

//------------------------------------------------------------------------------------------------------

} // namespace meshgen
} // namespace atlas

#endif // ReducedGridMeshGenerator_h
