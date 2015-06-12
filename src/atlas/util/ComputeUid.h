/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_ComputeUid_h
#define atlas_util_ComputeUid_h

#include <cmath>
#include <sstream>

#include "atlas/atlas_config.h"
#include "atlas/Parameters.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/LonLatPoint.h"

namespace atlas {
namespace util {

struct ComputeUid
{
  ComputeUid() {}
  ComputeUid( const FunctionSpace& nodes ):
    funcspace(&nodes)
  {
    update();
  }

  gidx_t operator()( const double crd[] ) const
  {
    return LonLatPoint( crd ).uid();
  }

  gidx_t operator()( int node ) const
  {
    return LonLatPoint( lonlat[node] ).uid();
  }

  gidx_t operator()( const IndexView<int,1>& elem_nodes ) const
  {
    double centroid[2];
    centroid[LON] = 0.;
    centroid[LAT] = 0.;
    int nb_elem_nodes = elem_nodes.shape(0);
    for( int jnode=0; jnode<nb_elem_nodes; ++jnode )
    {
      centroid[LON] += lonlat( elem_nodes(jnode), LON );
      centroid[LAT] += lonlat( elem_nodes(jnode), LAT );
    }
    centroid[LON] /= static_cast<double>(nb_elem_nodes);
    centroid[LAT] /= static_cast<double>(nb_elem_nodes);
    return LonLatPoint( centroid[LON], centroid[LAT] ).uid32();
  }

  gidx_t operator()( double crds[], int npts ) const
  {
    double centroid[2];
    centroid[LON] = 0.;
    centroid[LAT] = 0.;
    for( int jnode=0; jnode<npts; ++jnode )
    {
      centroid[LON] += crds[jnode*2+LON];
      centroid[LAT] += crds[jnode*2+LAT];
    }
    centroid[LON] /= static_cast<double>(npts);
    centroid[LAT] /= static_cast<double>(npts);
    return LonLatPoint( centroid ).uid32();
  }


  void update()
  {
    lonlat = ArrayView<double,2> ( funcspace->field("lonlat") );
  }
private:
  const FunctionSpace* funcspace;
  ArrayView<double,2> lonlat;
};

} // namespace util
} // namespace atlas

#endif
