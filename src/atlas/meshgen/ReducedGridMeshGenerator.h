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

namespace grids { class ReducedGrid; }
namespace meshgen {

struct Region;

class ReducedGridMeshGenerator
{
public:
  ReducedGridMeshGenerator();

  void generate( const grids::ReducedGrid&, Mesh& );

  Mesh* generate( const grids::ReducedGrid& );

  Mesh* operator()( const grids::ReducedGrid& );

private:
  void generate_region( const grids::ReducedGrid&, const std::vector<int>& parts, int mypart, Region& region );

  void generate_mesh( const grids::ReducedGrid&,const std::vector<int>& parts, const Region& region, Mesh& m );

  void generate_global_element_numbering( Mesh& mesh );

public:

  Metadata options;
  double max_angle_;
  bool triangulate_quads_;
};

} // namespace meshgen
} // namespace atlas

#endif // ReducedGridMeshGenerator_h
