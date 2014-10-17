/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef RGGMeshGenerator_h
#define RGGMeshGenerator_h

#include <vector>
#include <eckit/config/Configurable.h>
#include <eckit/config/Resource.h>
#include "atlas/Metadata.h"

namespace atlas {
  class Mesh;
namespace meshgen {

struct Region;
class RGG;

class RGGMeshGenerator
{
public:
  RGGMeshGenerator();
  Mesh* generate(const RGG& rgg);
  Mesh* operator()(const RGG& rgg){ return generate(rgg); }
private:
  void generate_region(const RGG& rgg, const std::vector<int>& parts, int mypart, Region& region);
  Mesh* generate_mesh(const RGG& rgg,const std::vector<int>& parts, const Region& region);
  void generate_global_element_numbering( Mesh& mesh );
  //std::vector<int> partition(const RGG& rgg) const;
public:
  Metadata options;
  double max_angle_;
	bool triangulate_quads_;
};

} // namespace meshgen
} // namespace atlas

#endif // RGGMeshGenerator_h
