/*
 * (C) Copyright 1996-2016 ECMWF.
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

namespace eckit { class Parametrisation; }

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
  ReducedGridMeshGenerator(const eckit::Parametrisation&);

  virtual void generate( const Grid&, const GridDistribution&, Mesh& ) const;
  virtual void generate( const Grid&, Mesh& ) const;
  using MeshGenerator::generate;

private:

  void configure_defaults();

  void generate_region( const grids::ReducedGrid&, const std::vector<int>& parts, int mypart, Region& region ) const;

  void generate_mesh_new( const grids::ReducedGrid&,const std::vector<int>& parts, const Region& region, Mesh& m ) const;

#if !DEPRECATE_OLD_FUNCTIONSPACE
  void generate_mesh_convert_to_old( const grids::ReducedGrid&,const std::vector<int>& parts, const Region& region, Mesh& m ) const;
#endif
  void generate_mesh( const grids::ReducedGrid&,const std::vector<int>& parts, const Region& region, Mesh& m ) const;

  void generate_global_element_numbering( Mesh& mesh ) const;

public:

  Metadata options;

};

//------------------------------------------------------------------------------------------------------

} // namespace meshgen
} // namespace atlas

#endif // ReducedGridMeshGenerator_h
