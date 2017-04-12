/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/mesh/actions/ExtendNodesGlobal.h"

#include "eckit/exception/Exceptions.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/util/CoordinateEnums.h"


namespace atlas {
namespace mesh {
namespace actions {

ExtendNodesGlobal::ExtendNodesGlobal( const std::string& gridname ) :
    gridname_(gridname) {
}


void ExtendNodesGlobal::operator()(const atlas::grid::Grid& grid, Mesh& mesh) const {

    if (grid.domain().global()) return; // don't add virtual points to global domains

    grid::Grid O16( "O16" );

    // virtual points
    std::vector<PointXY> extended_pts;
    extended_pts.reserve( grid.size() );

    // loop over the point and keep the ones that *don't* fall in the domain

    for( const PointLonLat& lonlat : O16.lonlat() ) {
      PointXY xy = grid.projection().xy(lonlat);
      if( not grid.domain().contains( xy ) ) {
        extended_pts.push_back(xy);
      }
    }

    mesh::Nodes& nodes = mesh.nodes();

    const size_t nb_real_pts = nodes.size();
    const size_t nb_extension_pts = extended_pts.size();

    size_t new_size = nodes.size() + extended_pts.size();

    nodes.resize(new_size); // resizes the fields

    const size_t nb_total_pts = nodes.size();

    ASSERT( nb_total_pts == nb_real_pts + nb_extension_pts );

    nodes.metadata().set<size_t>("NbRealPts",nb_real_pts);
    nodes.metadata().set<size_t>("NbVirtualPts",nb_extension_pts);

    array::ArrayView<double,2> xyz    = array::make_view<double,2>( nodes.field("xyz") );
    array::ArrayView<double,2> xy     = array::make_view<double,2>( nodes.xy() );
    array::ArrayView<double,2> lonlat = array::make_view<double,2>( nodes.lonlat() );
    array::ArrayView<gidx_t,1> gidx   = array::make_view<gidx_t,1>( nodes.global_index() );

    for(size_t i = 0; i < nb_extension_pts; ++i) {
        const size_t n = nb_real_pts + i;
        PointLonLat pLL  = grid.projection().lonlat(extended_pts[i]);
        PointXYZ    pXYZ = lonlat_to_geocentric(pLL);
        xyz(n,XX) = pXYZ.x();
        xyz(n,YY) = pXYZ.y();
        xyz(n,ZZ) = pXYZ.z();
        xy(n,XX) = extended_pts[i].x();
        xy(n,YY) = extended_pts[i].y();
        lonlat(n,LON) = pLL.lon();
        lonlat(n,LAT) = pLL.lat();
        gidx(n) = n+1;
    }

}


} // namespace actions
} // namespace mesh
} // namespace atlas
