/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <algorithm>

#define BOOST_TEST_MODULE TestHaloEastWest
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"

#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"
#include "atlas/meshgen/RGG.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Metadata.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/actions/BuildHalo.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/mesh/Parameters.hpp"

using namespace atlas;
using namespace atlas::meshgen;

#define DISABLE if(0)
#define ENABLE if(1)

namespace atlas {
namespace test {

class DebugMesh:   public RGG { public: DebugMesh(); };
DebugMesh::DebugMesh()
{
  int nlat=6;
  int lon[] = {
    12,
    12,
    12,
    12,
    12,
    12
  };
  /*
  First prediction of colatitudes
  */
  std::vector<double> colat(nlat);
  double z;
  for( int i=0; i<nlat; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*nlat+2.);
    colat[i] = z+1./(tan(z)*(8.*(2.*nlat)*(2.*nlat)));
  }
  /*
  Fill in final structures
  */
  lat_.resize(2*nlat);
  lon_.resize(2*nlat);
  std::copy( lon, lon+nlat, lon_.begin() );
  std::reverse_copy( lon, lon+nlat, lon_.begin()+nlat );
  std::copy( colat.begin(), colat.begin()+nlat, lat_.begin() );
  std::reverse_copy( colat.begin(), colat.begin()+nlat, lat_.begin()+nlat );
  for (int i=0; i<nlat; ++i)
    lat_[i]=M_PI/2.-lat_[i];
  for (int i=nlat; i<2*nlat; ++i)
    lat_[i]=-M_PI/2.+lat_[i];
}

} // end namespace test
} // end namespace atlas

BOOST_AUTO_TEST_CASE( init ) { MPL::init(); }

BOOST_AUTO_TEST_CASE( test_halo_2parts )
{
  RGGMeshGenerator generate;

  generate.options.set("nb_parts",MPL::size());
  generate.options.set("part",MPL::rank());

  Mesh* m = generate( T159() );
//  Mesh* m = generate( test::DebugMesh() );

  actions::build_parallel_fields(*m);
  actions::make_periodic(*m);
  actions::build_halo(*m,2);

  std::stringstream filename; filename << "T63_halo_EW_p"<<MPL::rank()<<".msh";
  Gmsh::write(*m,filename.str());
  delete m;
}

BOOST_AUTO_TEST_CASE( finalize ) { MPL::finalize(); }
