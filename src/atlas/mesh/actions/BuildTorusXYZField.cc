/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include "eckit/geometry/Point3.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildTorusXYZField.h"
#include "atlas/field/Field.h"
#include "atlas/array/ArrayView.h"
#include "atlas/grid/domain/RectangularDomain.h"
#include "atlas/grid/domain/GlobalDomain.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

BuildTorusXYZField::BuildTorusXYZField(const std::string& name)
    : name_(name)
{
}

field::Field& BuildTorusXYZField::operator()(Mesh& mesh, const atlas::grid::domain::Domain * dom, double r0, double r1, int nx, int ny) const
{
  return operator()(mesh.nodes(),dom,r0,r1,nx,ny);
}

field::Field& BuildTorusXYZField::operator()(mesh::Nodes& nodes, const atlas::grid::domain::Domain * dom, double r0, double r1, int nx, int ny) const
{
  // fill xyz with torus coordinates. r0 and r1 are large and small radii, respectively.

  std::vector<double> bbox;

  // check if the domain is rectangular
  const atlas::grid::domain::RectangularDomain* rdom=dynamic_cast<const atlas::grid::domain::RectangularDomain*>(dom);
  const atlas::grid::domain::GlobalDomain* gdom=dynamic_cast<const atlas::grid::domain::GlobalDomain*>(dom);
  if( rdom ) {
    bbox = rdom->bbox();
  }
  else if( gdom ) {
    bbox.resize(4);
    bbox[0]= 0.;
    bbox[1]= 360.;
    bbox[2]=-90.;
    bbox[3]= 90.;
  }
  if( bbox.empty() ) {
    throw eckit::BadCast("Torus can only be built from rectangular or global domain");
  }

  if( !nodes.has_field(name_) )
  {
    size_t npts = nodes.size();
    array::ArrayView<double,2> lonlat( nodes.lonlat() );
    array::ArrayView<double,2> xyz   ( nodes.add( field::Field::create<double>(name_,array::make_shape(npts,3) ) ) );

    const double pi=M_PI;
    const double c1 = 2.*pi/double(nx)*(nx-1)/(bbox[1]-bbox[0]);
    const double c2 = 2.*pi/double(ny)*(ny-1)/(bbox[3]-bbox[2]);
    for( size_t n=0; n<npts; ++n )
    {
      double *xx=xyz[n].data();
      double *ll=lonlat[n].data();

      double lon=-pi+c1*(ll[0]-bbox[0]);
      double lat=-pi+c2*(ll[1]-bbox[2]);

      xx[0]=std::cos(lon)*(r0+r1*std::cos(lat));
      xx[1]=std::sin(lon)*(r0+r1*std::cos(lat));
      xx[2]=r1*std::sin(lat);
    }
  }
  return nodes.field(name_);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas
