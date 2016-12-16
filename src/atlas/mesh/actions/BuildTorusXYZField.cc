/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/geometry/Point3.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildTorusXYZField.h"
#include "atlas/field/Field.h"
#include "atlas/array/ArrayView.h"

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


	// check if the domain is rectangular
	const atlas::grid::domain::RectangularDomain * rdom=dynamic_cast<const atlas::grid::domain::RectangularDomain *>(dom);
  if( !rdom )
    throw eckit::BadCast("Torus can only be built from rectangular domain",Here());
	
  if( !nodes.has_field(name_) )
  {
    size_t npts = nodes.size();
    array::ArrayView<double,2> lonlat( nodes.lonlat() );
    array::ArrayView<double,2> xyz   ( nodes.add( field::Field::create<double>(name_,array::make_shape(npts,3) ) ) );
    std::vector<double> bbox=rdom->bbox();
   	double pi=acos(-1.);
    for( size_t n=0; n<npts; ++n )
    {
    	double *xx=xyz[n].data();
    	double *ll=lonlat[n].data();
    	
    	double lon, lat;
    	lon=-pi+2*pi*(nx-1)*(ll[0]-bbox[0])/(bbox[1]-bbox[0])/nx;
    	lat=-pi+2*pi*(ny-1)*(ll[1]-bbox[2])/(bbox[3]-bbox[2])/ny;
    	
    	xx[0]=cos(lon)*(r0+r1*cos(lat));
    	xx[1]=sin(lon)*(r0+r1*cos(lat));
    	xx[2]=r1*sin(lat);
			
    }
  }
  return nodes.field(name_);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas
