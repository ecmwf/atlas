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
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/field/Field.h"
#include "atlas/array/ArrayView.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

BuildXYZField::BuildXYZField(const std::string& name)
    : name_(name)
{
}

field::Field& BuildXYZField::operator()(Mesh& mesh) const
{
  return operator()(mesh.nodes());
}

field::Field& BuildXYZField::operator()(mesh::Nodes& nodes) const
{
  if( !nodes.has_field(name_) )
  {
    size_t npts = nodes.size();
    array::ArrayView<double,2> lonlat( nodes.lonlat() );
    array::ArrayView<double,2> xyz   ( nodes.add( field::Field::create<double>(name_,array::make_shape(npts,3) ) ) );
    for( size_t n=0; n<npts; ++n )
    {
      eckit::geometry::lonlat_to_3d(lonlat[n].data(),xyz[n].data());
    }
  }
  return nodes.field(name_);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas
