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

#include "atlas/array.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/trans/Trans.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

#include "tests/AtlasTestEnvironment.h"

#if ATLAS_HAVE_TRANS
#include "atlas/trans/ifs/TransIFS.h"
#include "transi/trans.h"
#endif

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

struct AtlasTransEnvironment : public AtlasTestEnvironment {
    AtlasTransEnvironment( int argc, char* argv[] ) : AtlasTestEnvironment( argc, argv ) {
        if ( mpi::comm().size() == 1 ) {
            trans_use_mpi( false );
        }
        trans_init();
    }

    ~AtlasTransEnvironment() { trans_finalize(); }
};

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
static void rotated_flow_magnitude( StructuredGrid& grid, double var[], const double& beta ) {
    const double radius  = util::Earth::radius();
    const double USCAL   = 20.;
    const double pvel    = USCAL / radius;
    const double deg2rad = M_PI / 180.;

    idx_t n( 0 );
    for ( idx_t jlat = 0; jlat < grid.ny(); ++jlat ) {
        for ( idx_t jlon = 0; jlon < grid.nx( jlat ); ++jlon ) {
            const double x = grid.x( jlon, jlat ) * deg2rad;
            const double y = grid.y( jlat ) * deg2rad;
            const double Ux =
                pvel * ( std::cos( beta ) + std::tan( y ) * std::cos( x ) * std::sin( beta ) ) * radius * std::cos( y );
            const double Uy = -pvel * std::sin( x ) * std::sin( beta ) * radius;
            var[n]          = std::sqrt( Ux * Ux + Uy * Uy );
            ++n;
        }
    }
}

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
void rotated_flow_magnitude( const functionspace::NodeColumns& fs, Field& field, const double& beta ) {
    const double radius  = util::Earth::radius();
    const double USCAL   = 20.;
    const double pvel    = USCAL / radius;
    const double deg2rad = M_PI / 180.;

    array::ArrayView<double, 2> lonlat_deg = array::make_view<double, 2>( fs.nodes().lonlat() );
    array::ArrayView<double, 1> var        = array::make_view<double, 1>( field );

    size_t nnodes = fs.nodes().size();
    for ( size_t jnode = 0; jnode < nnodes; ++jnode ) {
        double x = lonlat_deg( jnode, (size_t)LON ) * deg2rad;
        double y = lonlat_deg( jnode, (size_t)LAT ) * deg2rad;
        double Ux =
            pvel * ( std::cos( beta ) + std::tan( y ) * std::cos( x ) * std::sin( beta ) ) * radius * std::cos( y );
        double Uy    = -pvel * std::sin( x ) * std::sin( beta ) * radius;
        var( jnode ) = std::sqrt( Ux * Ux + Uy * Uy );
    }
}

//-----------------------------------------------------------------------------

CASE( "test_invtrans_ifsStyle" ) {
    std::string grid_uid( "O80" );
    StructuredGrid g( grid_uid );
    long N = g.ny() / 2;
    trans::TransIFS trans( g, 2 * N - 1 );
    Log::info() << "Trans initialized" << std::endl;
    std::vector<double> rspecg;
    int nfld = 1;

    std::vector<double> init_gpg( trans.grid().size() );
    std::vector<double> init_gp( trans.trans()->ngptot );
    std::vector<double> init_sp( trans.trans()->nspec2 );
    std::vector<int> nfrom( nfld, 1 );
    if ( mpi::comm().rank() == 0 ) {
        double beta = M_PI * 0.5;
        rotated_flow_magnitude( g, init_gpg.data(), beta );
    }
    trans.distgrid( nfld, nfrom.data(), init_gpg.data(), init_gp.data() );
    trans.dirtrans( nfld, init_gp.data(), init_sp.data() );

    std::vector<double> rgp( 3 * nfld * trans.trans()->ngptot );
    double *no_vorticity( nullptr ), *no_divergence( nullptr );
    int nb_vordiv( 0 );
    int nb_scalar( nfld );
    trans.invtrans( nb_scalar, init_sp.data(), nb_vordiv, no_vorticity, no_divergence, rgp.data(),
                    option::scalar_derivatives( true ) );

    std::vector<int> nto( nfld, 1 );
    std::vector<double> rgpg( 3 * nfld * trans.grid().size() );

    trans.gathgrid( nfld, nto.data(), rgp.data(), rgpg.data() );

    // Output
    {
        Mesh mesh = StructuredMeshGenerator().generate( g );
        functionspace::StructuredColumns gp( g );
        output::Gmsh gmsh( grid_uid + "-grid.msh" );
        Field scalar( "scalar", rgp.data(), array::make_shape( gp.size() ) );
        Field scalar_dNS( "scalar_dNS", rgp.data() + nfld * gp.size(), array::make_shape( gp.size() ) );
        Field scalar_dEW( "scalar_dEW", rgp.data() + 2 * nfld * gp.size(), array::make_shape( gp.size() ) );
        gmsh.write( mesh );
        gmsh.write( scalar, gp );
        gmsh.write( scalar_dEW, gp );
        gmsh.write( scalar_dNS, gp );
    }
}

CASE( "test_invtrans_grad" ) {
    std::string grid_uid( "O48" );
    StructuredGrid g( grid_uid );
    Mesh mesh = StructuredMeshGenerator().generate( g );
    idx_t N   = g.ny() / 2;
    trans::Trans trans( g, 2 * N - 1 );
    functionspace::NodeColumns gp( mesh );
    functionspace::Spectral sp( trans );

    Field scalar_sp = sp.createField<double>( option::name( "scalar_sp" ) );
    Field scalar    = gp.createField<double>( option::name( "scalar" ) );
    Field grad      = gp.createField<double>( option::name( "grad" ) | option::variables( 2 ) );

    // Initial condition
    double beta = M_PI * 0.5;
    rotated_flow_magnitude( gp, scalar, beta );

    // Transform to spectral
    trans.dirtrans( scalar, scalar_sp );

    // Inverse transform for gradient
    trans.invtrans_grad( scalar_sp, grad );

    gp.haloExchange( grad );

    // Output
    {
        Mesh mesh = StructuredMeshGenerator().generate( g );
        functionspace::StructuredColumns gp( g );
        output::Gmsh gmsh( grid_uid + "-nodes.msh" );
        gmsh.write( mesh );
        gmsh.write( scalar, gp );
        gmsh.write( grad, gp );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run<atlas::test::AtlasTransEnvironment>( argc, argv );
}
