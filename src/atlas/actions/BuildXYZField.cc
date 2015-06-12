
#include "atlas/actions/BuildXYZField.h"
#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/util/ArrayView.h"
#include "eckit/geometry/Point3.h"

namespace atlas {
namespace actions {

Field& build_xyz_field( Mesh& mesh,  const std::string& name)
{
  return build_xyz_field( mesh.function_space("nodes"), name );
}

Field& build_xyz_field( FunctionSpace& nodes, const std::string& name )
{
  if( !nodes.has_field(name) )
  {
    ArrayView<double,2> lonlat( nodes.field("lonlat") );
    ArrayView<double,2> xyz   ( nodes.create_field<double>(name,3) );
    size_t npts = nodes.shape(0);
    for( size_t n=0; n<npts; ++n )
    {
      eckit::geometry::lonlat_to_3d(lonlat[n].data(),xyz[n].data());
    }
  }
  return nodes.field(name);
}

}
}
