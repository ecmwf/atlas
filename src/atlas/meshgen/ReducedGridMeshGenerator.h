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

namespace atlas {
class Mesh;
class GridDistribution;

namespace grids { class ReducedGrid; }
namespace meshgen {

struct Region;

//------------------------------------------------------------------------------------------------------

class ReducedGridMeshGenerator {

public:

  ReducedGridMeshGenerator();

  void generate( const grids::ReducedGrid&, const GridDistribution&, Mesh& );
  void generate( const grids::ReducedGrid&, GridDistribution*, Mesh& );
  void generate( const grids::ReducedGrid&, Mesh& );

  Mesh* generate( const grids::ReducedGrid&, const GridDistribution& );
  Mesh* generate( const grids::ReducedGrid&, GridDistribution* );
  Mesh* generate( const grids::ReducedGrid& );

  Mesh* operator()( const grids::ReducedGrid&, const GridDistribution& );
  Mesh* operator()( const grids::ReducedGrid&, GridDistribution* );
  Mesh* operator()( const grids::ReducedGrid& );


  void set_three_dimensional( bool );
  void set_patch_pole( bool );
  void set_include_pole( bool );

private:

  void generate_region( const grids::ReducedGrid&, const std::vector<int>& parts, int mypart, Region& region );

  void generate_mesh( const grids::ReducedGrid&,const std::vector<int>& parts, const Region& region, Mesh& m );

  void generate_global_element_numbering( Mesh& mesh );

public:

  Metadata options;

};

//------------------------------------------------------------------------------------------------------

} // namespace meshgen
} // namespace atlas

#endif // ReducedGridMeshGenerator_h
