/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <iostream>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/numerics/Nabla.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::numerics;
using namespace atlas::meshgenerator;
using namespace atlas::grid;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

double dual_volume( const Mesh& mesh ) {
    const mesh::Nodes& nodes = mesh.nodes();
    int nb_nodes             = nodes.size();
    auto dual_volumes        = array::make_view<double, 1>( nodes.field( "dual_volumes" ) );
    auto is_ghost            = array::make_view<int, 1>( nodes.ghost() );
    double area              = 0;
    for ( int node = 0; node < nb_nodes; ++node ) {
        if ( !is_ghost( node ) ) {
            area += dual_volumes( node );
        }
    }

    mpi::comm().allReduceInPlace( area, eckit::mpi::sum() );

    return area;
}

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
void rotated_flow( const fvm::Method& fvm, Field& field, const double& beta ) {
    const double radius  = fvm.radius();
    const double USCAL   = 20.;
    const double pvel    = USCAL / radius;
    const double deg2rad = M_PI / 180.;

    auto lonlat_deg = array::make_view<double, 2>( fvm.mesh().nodes().lonlat() );
    auto var        = array::make_view<double, 3>( field );

    idx_t nnodes = fvm.mesh().nodes().size();
    for ( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        double x = lonlat_deg( jnode, LON ) * deg2rad;
        double y = lonlat_deg( jnode, LAT ) * deg2rad;
        double Ux =
            pvel * ( std::cos( beta ) + std::tan( y ) * std::cos( x ) * std::sin( beta ) ) * radius * std::cos( y );
        double Uy = -pvel * std::sin( x ) * std::sin( beta ) * radius;
        for ( idx_t jlev = 0; jlev < field.levels(); ++jlev ) {
            var( jnode, jlev, LON ) = Ux;
            var( jnode, jlev, LAT ) = Uy;
        }
    }
}

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
void rotated_flow_magnitude( const fvm::Method& fvm, Field& field, const double& beta ) {
    const double radius  = fvm.radius();
    const double USCAL   = 20.;
    const double pvel    = USCAL / radius;
    const double deg2rad = M_PI / 180.;

    auto lonlat_deg = array::make_view<double, 2>( fvm.mesh().nodes().lonlat() );
    auto var        = array::make_view<double, 2>( field );

    idx_t nnodes = fvm.mesh().nodes().size();
    for ( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        double x = lonlat_deg( jnode, LON ) * deg2rad;
        double y = lonlat_deg( jnode, LAT ) * deg2rad;
        double Ux =
            pvel * ( std::cos( beta ) + std::tan( y ) * std::cos( x ) * std::sin( beta ) ) * radius * std::cos( y );
        double Uy = -pvel * std::sin( x ) * std::sin( beta ) * radius;
        for ( idx_t jlev = 0; jlev < field.levels(); ++jlev ) {
            var( jnode, jlev ) = std::sqrt( Ux * Ux + Uy * Uy );
        }
    }
}

static std::string griduid() {
    return "Slat20";
}

//-----------------------------------------------------------------------------

CASE( "test_factory" ) {
    EXPECT( NablaFactory::has( "fvm" ) );
}

CASE( "test_build" ) {
    Log::info() << "test_build" << std::endl;
    MeshGenerator meshgenerator( "structured" );
    Mesh mesh      = meshgenerator.generate( Grid( "O16" ) );
    const double R = util::Earth::radius();
    fvm::Method fvm( mesh, util::Config( "radius", R ) );
    Nabla nabla( fvm );

    double spherical_area = 360. * 180.;
    EXPECT( eckit::types::is_approximately_equal( dual_volume( mesh ), spherical_area, 5.0 ) );
}

CASE( "test_grad" ) {
    Log::info() << "test_grad" << std::endl;
    idx_t nlev  = 1;
    auto radius = option::radius( "Earth" );
    Grid grid( griduid() );
    MeshGenerator meshgenerator( "structured" );
    Mesh mesh = meshgenerator.generate( grid, Distribution( grid, Partitioner( "equal_regions" ) ) );
    fvm::Method fvm( mesh, radius | option::levels( nlev ) );
    Nabla nabla( fvm );

    idx_t nnodes = mesh.nodes().size();

    FieldSet fields;
    fields.add( fvm.node_columns().createField<double>( option::name( "scalar" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "rscalar" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "grad" ) | option::variables( 2 ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "rgrad" ) | option::variables( 2 ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "xder" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "yder" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "rxder" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "ryder" ) ) );

    EXPECT( fields["scalar"].rank() == 2 );
    EXPECT( fields["grad"].rank() == 3 );
    EXPECT( fields["scalar"].levels() == nlev );
    EXPECT( fields["grad"].levels() == nlev );
    //  fields.add( fvm.createField<double>("exact_yder",nlev) );

    //  const double deg2rad = M_PI/180.;
    //  array::ArrayView<double,2> var( fields["scalar"] );
    ////  array::ArrayView<double,2> exact_yder( fields["exact_yder"] );
    //  for( size_t jnode=0; jnode< nnodes ; ++jnode )
    //  {
    //    const double y  = lonlat(jnode,LAT) * deg2rad ;

    //    for(size_t jlev = 0; jlev < nlev; ++jlev) {
    //      var(jnode,jlev)        = std::sin(4.*y);
    ////      exact_yder(jnode,jlev) = 4.*std::cos(4.*y)/radius;
    //    }
    //  }

    rotated_flow_magnitude( fvm, fields["scalar"], 0. );
    rotated_flow_magnitude( fvm, fields["rscalar"], M_PI_2 * 0.75 );

    nabla.gradient( fields["scalar"], fields["grad"] );
    nabla.gradient( fields["rscalar"], fields["rgrad"] );
    auto xder        = array::make_view<double, 2>( fields["xder"] );
    auto yder        = array::make_view<double, 2>( fields["yder"] );
    auto rxder       = array::make_view<double, 2>( fields["rxder"] );
    auto ryder       = array::make_view<double, 2>( fields["ryder"] );
    const auto grad  = array::make_view<double, 3>( fields["grad"] );
    const auto rgrad = array::make_view<double, 3>( fields["rgrad"] );
    for ( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        for ( idx_t jlev = 0; jlev < nlev; ++jlev ) {
            xder( jnode, jlev )  = grad( jnode, jlev, LON );
            yder( jnode, jlev )  = grad( jnode, jlev, LAT );
            rxder( jnode, jlev ) = rgrad( jnode, jlev, LON );
            ryder( jnode, jlev ) = rgrad( jnode, jlev, LAT );
        }
    }

    // output to gmsh
    {
        fvm.node_columns().haloExchange( fields );
        output::Gmsh( grid.name() + ".msh" ).write( mesh );
        output::Gmsh gmsh_fields( grid.name() + "_fields.msh" );
        gmsh_fields.write( fields["scalar"] );
        gmsh_fields.write( fields["xder"] );
        gmsh_fields.write( fields["yder"] );
        gmsh_fields.write( fields["rscalar"] );
        gmsh_fields.write( fields["rxder"] );
        gmsh_fields.write( fields["ryder"] );
    }
}

CASE( "test_div" ) {
    Log::info() << "test_div" << std::endl;
    size_t nlev         = 1;
    const double radius = util::Earth::radius();
    //  const double radius = 1.;
    Grid grid( griduid() );
    MeshGenerator meshgenerator( "structured" );
    Mesh mesh = meshgenerator.generate( grid, Distribution( grid, Partitioner( "equal_regions" ) ) );
    fvm::Method fvm( mesh, util::Config( "radius", radius ) | option::levels( nlev ) );
    Nabla nabla( fvm );

    FieldSet fields;
    fields.add( fvm.node_columns().createField<double>( option::name( "wind" ) | option::variables( 2 ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "div" ) ) );

    rotated_flow( fvm, fields["wind"], M_PI_2 * 0.75 );

    nabla.divergence( fields["wind"], fields["div"] );

    // output to gmsh
    {
        fvm.node_columns().haloExchange( fields );
        output::Gmsh gmsh( grid.name() + "_fields.msh", "a" );
        gmsh.write( fields["wind"] );
        gmsh.write( fields["div"] );
    }
}

CASE( "test_curl" ) {
    Log::info() << "test_curl" << std::endl;
    size_t nlev         = 1;
    const double radius = util::Earth::radius();
    //  const double radius = 1.;
    Grid grid( griduid() );
    MeshGenerator meshgenerator( "structured" );
    Mesh mesh = meshgenerator.generate( grid, Distribution( grid, Partitioner( "equal_regions" ) ) );
    fvm::Method fvm( mesh, util::Config( "radius", radius ) | option::levels( nlev ) );
    Nabla nabla( fvm );

    FieldSet fields;
    fields.add( fvm.node_columns().createField<double>( option::name( "wind" ) | option::variables( 2 ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "vor" ) ) );

    rotated_flow( fvm, fields["wind"], M_PI_2 * 0.75 );

    nabla.curl( fields["wind"], fields["vor"] );

    fields.add( fvm.node_columns().createField<double>( option::name( "windgrad" ) | option::variables( 2 * 2 ) ) );
    nabla.gradient( fields["wind"], fields["windgrad"] );

    fields.add( fvm.node_columns().createField<double>( option::name( "windX" ) | option::levels( false ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "windY" ) | option::levels( false ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "windXgradX" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "windXgradY" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "windYgradX" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "windYgradY" ) ) );
    auto wind     = array::make_view<double, 3>( fields["wind"] );
    auto windgrad = array::make_view<double, 3>( fields["windgrad"] );

    auto windX      = array::make_view<double, 1>( fields["windX"] );
    auto windY      = array::make_view<double, 1>( fields["windY"] );
    auto windXgradX = array::make_view<double, 2>( fields["windXgradX"] );
    auto windXgradY = array::make_view<double, 2>( fields["windXgradY"] );
    auto windYgradX = array::make_view<double, 2>( fields["windYgradX"] );
    auto windYgradY = array::make_view<double, 2>( fields["windYgradY"] );
    for ( idx_t j = 0; j < windX.size(); ++j ) {
        static const idx_t lev0 = 0;
        static const idx_t XdX  = XX * 2 + XX;
        static const idx_t XdY  = XX * 2 + YY;
        static const idx_t YdX  = YY * 2 + XX;
        static const idx_t YdY  = YY * 2 + YY;
        windX( j )              = wind( j, lev0, XX );
        windY( j )              = wind( j, lev0, YY );
        windXgradX( j, lev0 )   = windgrad( j, lev0, XdX );
        windXgradY( j, lev0 )   = windgrad( j, lev0, XdY );
        windYgradX( j, lev0 )   = windgrad( j, lev0, YdX );
        windYgradY( j, lev0 )   = windgrad( j, lev0, YdY );
    }

    // output to gmsh
    {
        fvm.node_columns().haloExchange( fields );
        output::Gmsh gmsh( grid.name() + "_fields.msh", "a" );
        gmsh.write( fields["vor"] );
        gmsh.write( fields["windX"] );
        gmsh.write( fields["windXgradX"] );
        gmsh.write( fields["windXgradY"] );
        gmsh.write( fields["windY"] );
        gmsh.write( fields["windYgradX"] );
        gmsh.write( fields["windYgradY"] );
        gmsh.write( fields["windgrad"] );
    }
}

CASE( "test_lapl" ) {
    Log::info() << "test_lapl" << std::endl;
    size_t nlev         = 1;
    const double radius = util::Earth::radius();
    //  const double radius = 1.;
    Grid grid( griduid() );
    MeshGenerator meshgenerator( "structured" );
    Mesh mesh = meshgenerator.generate( grid, Distribution( grid, Partitioner( "equal_regions" ) ) );
    fvm::Method fvm( mesh, util::Config( "radius", radius ) | option::levels( nlev ) );
    Nabla nabla( fvm );

    FieldSet fields;
    fields.add( fvm.node_columns().createField<double>( option::name( "scal" ) ) );
    fields.add( fvm.node_columns().createField<double>( option::name( "lapl" ) ) );

    rotated_flow_magnitude( fvm, fields["scal"], M_PI_2 * 0.75 );

    nabla.laplacian( fields["scal"], fields["lapl"] );

    // output to gmsh
    {
        fvm.node_columns().haloExchange( fields );
        output::Gmsh gmsh( grid.name() + "_fields.msh", "a" );
        gmsh.write( fields["lapl"] );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
