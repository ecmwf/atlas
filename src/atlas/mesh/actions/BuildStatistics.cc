/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>

#include "eckit/filesystem/PathName.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/library/config.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"
#include "atlas/util/UnitSphere.h"

namespace atlas {
namespace mesh {
namespace actions {

namespace {

void tri_quality( double& eta, double& rho, const PointLonLat& p1, const PointLonLat& p2, const PointLonLat& p3 ) {
    // see
    // http://www.gidhome.com/component/manual/referencemanual/preprocessing/mesh_menu/mesh_quality

    double l12 = util::Constants::radiansToDegrees() * util::Earth::centralAngle( p1, p2 );
    double l23 = util::Constants::radiansToDegrees() * util::Earth::centralAngle( p2, p3 );
    double l31 = util::Constants::radiansToDegrees() * util::Earth::centralAngle( p3, p1 );

    double s    = 0.5 * ( l12 + l23 + l31 );
    double area = std::sqrt( s * ( s - l12 ) * ( s - l23 ) * ( s - l31 ) );

    eta = ( 4 * area * std::sqrt( 3. ) ) / ( std::pow( l12, 2 ) + std::pow( l23, 2 ) + std::pow( l31, 2 ) );

    double min_length = std::min( std::min( l12, l23 ), l31 );
    double max_length = std::max( std::max( l12, l23 ), l31 );

    rho = min_length / max_length;
}

void quad_quality( double& eta, double& rho, const PointLonLat& p1, const PointLonLat& p2, const PointLonLat& p3,
                   const PointLonLat& p4 ) {
    // see http://geuz.org/gmsh/doc/preprints/gmsh_quad_preprint.pdf

    PointXYZ xyz[4];
    util::UnitSphere::convertSphericalToCartesian( p1, xyz[0] );
    util::UnitSphere::convertSphericalToCartesian( p2, xyz[1] );
    util::UnitSphere::convertSphericalToCartesian( p3, xyz[2] );
    util::UnitSphere::convertSphericalToCartesian( p4, xyz[3] );

    PointXYZ l2m1( PointXYZ::sub( xyz[1], xyz[0] ) );
    PointXYZ l3m2( PointXYZ::sub( xyz[2], xyz[1] ) );
    PointXYZ l4m3( PointXYZ::sub( xyz[3], xyz[2] ) );
    PointXYZ l1m4( PointXYZ::sub( xyz[0], xyz[3] ) );

    double norm_l2m1 = PointXYZ::norm( l2m1 );
    double norm_l3m2 = PointXYZ::norm( l3m2 );
    double norm_l4m3 = PointXYZ::norm( l4m3 );
    double norm_l1m4 = PointXYZ::norm( l1m4 );

    double dot_l4m1_l2m1 = -PointXYZ::dot( l1m4, l2m1 );
    double dot_l1m2_l3m2 = -PointXYZ::dot( l2m1, l3m2 );
    double dot_l2m3_l4m3 = -PointXYZ::dot( l3m2, l4m3 );
    double dot_l3m4_l1m4 = -PointXYZ::dot( l4m3, l1m4 );

    // Angles at each quad corner
    double a1 = std::acos( dot_l4m1_l2m1 / ( norm_l1m4 * norm_l2m1 ) );
    double a2 = std::acos( dot_l1m2_l3m2 / ( norm_l2m1 * norm_l3m2 ) );
    double a3 = std::acos( dot_l2m3_l4m3 / ( norm_l3m2 * norm_l4m3 ) );
    double a4 = std::acos( dot_l3m4_l1m4 / ( norm_l4m3 * norm_l1m4 ) );

    double max_inner =
        std::max( std::max( std::max( std::abs( M_PI_2 - a1 ), std::abs( M_PI_2 - a2 ) ), std::abs( M_PI_2 - a3 ) ),
                  std::abs( M_PI_2 - a4 ) );

    eta = std::max( 1. - M_2_PI * max_inner, 0. );

    double l12 = util::Earth::centralAngle( p1, p2 );
    double l23 = util::Earth::centralAngle( p2, p3 );
    double l34 = util::Earth::centralAngle( p3, p4 );
    double l41 = util::Earth::centralAngle( p4, p1 );

    double min_length = std::min( std::min( std::min( l12, l23 ), l34 ), l41 );
    double max_length = std::max( std::max( std::max( l12, l23 ), l34 ), l41 );

    rho = min_length / max_length;
}

}  // namespace

void build_statistics( Mesh& mesh ) {
    mesh::Nodes& nodes                 = mesh.nodes();
    array::ArrayView<double, 2> lonlat = array::make_view<double, 2>( nodes.lonlat() );

    if ( mesh.edges().size() ) {
        if ( !mesh.edges().has_field( "arc_length" ) ) {
            mesh.edges().add(
                Field( "arc_length", array::make_datatype<double>(), array::make_shape( mesh.edges().size() ) ) );
        }
        array::ArrayView<double, 1> dist = array::make_view<double, 1>( mesh.edges().field( "arc_length" ) );
        const mesh::HybridElements::Connectivity& edge_nodes = mesh.edges().node_connectivity();

        const int nb_edges = mesh.edges().size();
        for ( int jedge = 0; jedge < nb_edges; ++jedge ) {
            int ip1 = edge_nodes( jedge, 0 );
            int ip2 = edge_nodes( jedge, 1 );
            PointLonLat p1( lonlat( ip1, LON ), lonlat( ip1, LAT ) );
            PointLonLat p2( lonlat( ip2, LON ), lonlat( ip2, LAT ) );
            dist( jedge ) = util::Earth::distance( p1, p2 ) * 1e-3;
        }
    }

    std::ofstream ofs;
    eckit::PathName stats_path( "stats.txt" );
    int idt = 10;
    if ( mpi::size() == 1 ) {
        ofs.open( stats_path.localPath(), std::ofstream::out );
        ofs << "# STATISTICS rho (min_length/max_length), eta (quality) \n";
        ofs << std::setw( idt ) << "# rho";
        ofs << std::setw( idt ) << "eta";
        ofs << "\n";
        ofs.close();
    }

    // Cell statistics
    {
        if ( mpi::size() == 1 ) {
            ofs.open( stats_path.localPath(), std::ofstream::app );
        }

        array::ArrayView<double, 1> rho = array::make_view<double, 1>( mesh.cells().add(
            Field( "stats_rho", array::make_datatype<double>(), array::make_shape( mesh.cells().size() ) ) ) );
        array::ArrayView<double, 1> eta = array::make_view<double, 1>( mesh.cells().add(
            Field( "stats_eta", array::make_datatype<double>(), array::make_shape( mesh.cells().size() ) ) ) );

        for ( idx_t jtype = 0; jtype < mesh.cells().nb_types(); ++jtype ) {
            const mesh::Elements& elements      = mesh.cells().elements( jtype );
            const BlockConnectivity& elem_nodes = elements.node_connectivity();
            const idx_t nb_elems                = elements.size();

            if ( elements.element_type().name() == "Triangle" ) {
                for ( idx_t jelem = 0; jelem < nb_elems; ++jelem ) {
                    idx_t ielem = elements.begin() + jelem;
                    idx_t ip1   = elem_nodes( jelem, 0 );
                    idx_t ip2   = elem_nodes( jelem, 1 );
                    idx_t ip3   = elem_nodes( jelem, 2 );
                    PointLonLat p1( lonlat( ip1, LON ), lonlat( ip1, LAT ) );
                    PointLonLat p2( lonlat( ip2, LON ), lonlat( ip2, LAT ) );
                    PointLonLat p3( lonlat( ip3, LON ), lonlat( ip3, LAT ) );

                    tri_quality( eta( ielem ), rho( ielem ), p1, p2, p3 );

                    if ( mpi::size() == 1 ) {
                        ofs << std::setw( idt ) << rho( ielem ) << std::setw( idt ) << eta( ielem ) << "\n";
                    }
                }
            }
            if ( elements.element_type().name() == "Quadrilateral" ) {
                for ( idx_t jelem = 0; jelem < nb_elems; ++jelem ) {
                    idx_t ielem = elements.begin() + jelem;
                    idx_t ip1   = elem_nodes( jelem, 0 );
                    idx_t ip2   = elem_nodes( jelem, 1 );
                    idx_t ip3   = elem_nodes( jelem, 2 );
                    idx_t ip4   = elem_nodes( jelem, 3 );

                    PointLonLat p1( lonlat( ip1, LON ), lonlat( ip1, LAT ) );
                    PointLonLat p2( lonlat( ip2, LON ), lonlat( ip2, LAT ) );
                    PointLonLat p3( lonlat( ip3, LON ), lonlat( ip3, LAT ) );
                    PointLonLat p4( lonlat( ip4, LON ), lonlat( ip4, LAT ) );

                    quad_quality( eta( ielem ), rho( ielem ), p1, p2, p3, p4 );

                    if ( mpi::size() == 1 ) {
                        ofs << std::setw( idt ) << rho( ielem ) << std::setw( idt ) << eta( ielem ) << "\n";
                    }
                }
            }
        }
        if ( mpi::size() == 1 ) {
            ofs.close();
        }
    }

    eckit::PathName dual_stats_path( "dual_stats.txt" );
    if ( mpi::size() == 1 ) {
        ofs.open( dual_stats_path.localPath(), std::ofstream::out );
        ofs << "# STATISTICS dual_area \n";
        ofs << std::setw( idt ) << "# area";
        ofs << "\n";
    }

    if ( nodes.has_field( "dual_volumes" ) ) {
        array::ArrayView<double, 1> dual_volumes   = array::make_view<double, 1>( nodes.field( "dual_volumes" ) );
        array::ArrayView<double, 1> dual_delta_sph = array::make_view<double, 1>( nodes.add(
            Field( "dual_delta_sph", array::make_datatype<double>(), array::make_shape( nodes.size(), 1 ) ) ) );

        for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
            const double lat        = util::Constants::degreesToRadians() * lonlat( jnode, LAT );
            const double hx         = util::Constants::degreesToRadians() * util::Earth::radius() * std::cos( lat );
            const double hy         = util::Constants::degreesToRadians() * util::Earth::radius();
            dual_delta_sph( jnode ) = std::sqrt( dual_volumes( jnode ) * hx * hy );
        }

        if ( mpi::size() == 1 ) {
            for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
                ofs << std::setw( idt ) << dual_delta_sph( jnode ) << "\n";
            }
        }
    }
    if ( mpi::size() == 1 ) {
        ofs.close();
    }
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_statistics( Mesh::Implementation* mesh ) {
    ATLAS_ASSERT( mesh != nullptr, "Cannot access uninitialised atlas_Mesh" );
    Mesh m( mesh );
    build_statistics( m );
}

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
