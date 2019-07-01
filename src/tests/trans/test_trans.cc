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

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/io/DataHandle.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/VorDivToUV.h"
#include "atlas/trans/detail/TransFactory.h"

#include "tests/AtlasTestEnvironment.h"

#if ATLAS_HAVE_TRANS
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/ifs/TransIFSNodeColumns.h"
#include "atlas/trans/ifs/TransIFSStructuredColumns.h"
#include "transi/trans.h"
#endif

using namespace eckit;
using atlas::grid::detail::partitioner::EqualRegionsPartitioner;
using atlas::grid::detail::partitioner::TransPartitioner;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

struct AtlasTransEnvironment : public AtlasTestEnvironment {
    AtlasTransEnvironment( int argc, char* argv[] ) : AtlasTestEnvironment( argc, argv ) {
        if ( mpi::comm().size() == 1 ) trans_use_mpi( false );
        trans_init();
    }

    ~AtlasTransEnvironment() { trans_finalize(); }
};

//-----------------------------------------------------------------------------

void read_rspecg( const trans::TransImpl& trans, std::vector<double>& rspecg, std::vector<int>& nfrom, int& nfld ) {
    Log::info() << "read_rspecg ...\n";
    nfld = 2;
    if ( mpi::comm().rank() == 0 ) {
        rspecg.resize( nfld * trans.nb_spectral_coefficients_global() );
        for ( size_t i = 0; i < trans.nb_spectral_coefficients_global(); ++i ) {
            rspecg[i * nfld + 0] = ( i == 0 ? 1. : 0. );  // scalar field 1
            rspecg[i * nfld + 1] = ( i == 0 ? 2. : 0. );  // scalar field 2
        }
    }
    nfrom.resize( nfld );
    for ( int jfld = 0; jfld < nfld; ++jfld )
        nfrom[jfld] = 1;

    Log::info() << "read_rspecg ... done" << std::endl;
}

//-----------------------------------------------------------------------------

void read_rspecg( Field spec ) {
    Log::info() << "read_rspecg ...\n";
    if ( mpi::comm().rank() == 0 ) {
        functionspace::Spectral funcspace = spec.functionspace();
        int nb_spectral_coefficients_global =
            functionspace::Spectral( spec.functionspace() ).nb_spectral_coefficients_global();
        auto view = array::make_view<double, 2>( spec );
        ATLAS_ASSERT( view.shape( 1 ) >= 2 );
        for ( int i = 0; i < nb_spectral_coefficients_global; ++i ) {
            view( i, 0 ) = ( i == 0 ? 1. : 0. );  // scalar field 1
            view( i, 1 ) = ( i == 0 ? 2. : 0. );  // scalar field 2
        }
    }
    Log::info() << "read_rspecg ... done" << std::endl;
}

//-----------------------------------------------------------------------------

CASE( "test_trans_distribution_matches_atlas" ) {
    EXPECT( grid::Partitioner::exists( "trans" ) );

    // Create grid and trans object
    Grid g( "N80" );

    EXPECT( StructuredGrid( g ).ny() == 160 );

    auto trans_partitioner = new TransPartitioner();
    grid::Partitioner partitioner( trans_partitioner );
    grid::Distribution distribution( g, partitioner );

    trans::TransIFS trans( g, 159 );
    ::Trans_t* t = trans;

    ATLAS_DEBUG_VAR( trans.truncation() );
    EXPECT( trans.truncation() == 159 );

    // -------------- do checks -------------- //
    EXPECT( t->nproc == int( mpi::comm().size() ) );
    EXPECT( t->myproc == int( mpi::comm().rank() + 1 ) );

    if ( mpi::comm().rank() == 0 )  // all tasks do the same, so only one needs to check
    {
        int max_nb_regions_EW( 0 );
        for ( int j = 0; j < trans_partitioner->nb_bands(); ++j )
            max_nb_regions_EW = std::max( max_nb_regions_EW, trans_partitioner->nb_regions( j ) );

        EXPECT( t->n_regions_NS == trans_partitioner->nb_bands() );
        EXPECT( t->n_regions_EW == max_nb_regions_EW );

        EXPECT( distribution.nb_partitions() == idx_t( mpi::comm().size() ) );
        EXPECT( idx_t( distribution.partition().size() ) == g.size() );

        std::vector<int> npts( distribution.nb_partitions(), 0 );

        for ( idx_t j = 0; j < g.size(); ++j )
            ++npts[distribution.partition( j )];

        EXPECT( t->ngptotg == g.size() );
        EXPECT( t->ngptot == npts[mpi::comm().rank()] );
        EXPECT( t->ngptotmx == *std::max_element( npts.begin(), npts.end() ) );

        // array::LocalView<int,1> n_regions ( trans.n_regions() ) ;
        for ( int j = 0; j < trans_partitioner->nb_bands(); ++j )
            EXPECT( t->n_regions[j] == trans_partitioner->nb_regions( j ) );
    }
}

CASE( "test_trans_options" ) {
    util::Config opts( option::fft( "FFTW" ) | option::split_latitudes( false ) | option::read_legendre( "readfile" ) );
    Log::info() << "trans_opts = " << opts << std::endl;
}

#ifdef TRANS_HAVE_IO
CASE( "test_write_read_cache" ) {
    Log::info() << "test_write_read_cache" << std::endl;
    using namespace trans;
    if ( mpi::comm().size() == 1 ) {
        // Create trans that will write file
        Trans trans_write_F24( Grid( "F24" ), 23,
                               option::write_legendre( "cached_legendre_coeffs-F24" ) | option::flt( false ) );
        Trans trans_write_N24( Grid( "N24" ), 23,
                               option::write_legendre( "cached_legendre_coeffs-N24" ) | option::flt( false ) );
        Trans trans_write_O24( Grid( "O24" ), 23,
                               option::write_legendre( "cached_legendre_coeffs-O24" ) | option::flt( false ) );

        // Create trans that will read from file
        Trans trans_read_F24( Grid( "F24" ), 23,
                              option::read_legendre( "cached_legendre_coeffs-F24" ) | option::flt( false ) );
        Trans trans_read_N24( Grid( "N24" ), 23,
                              option::read_legendre( "cached_legendre_coeffs-N24" ) | option::flt( false ) );
        Trans trans_read_O24( Grid( "O24" ), 23,
                              option::read_legendre( "cached_legendre_coeffs-O24" ) | option::flt( false ) );

        LegendreCache legendre_cache_F24( "cached_legendre_coeffs-F24" );
        LegendreCache legendre_cache_N24( "cached_legendre_coeffs-N24" );
        LegendreCache legendre_cache_O24( "cached_legendre_coeffs-O24" );

        Trans trans_cache_F24( legendre_cache_F24, Grid( "F24" ), 23, option::flt( false ) );
        Trans trans_cache_N24( legendre_cache_N24, Grid( "N24" ), 23, option::flt( false ) );
        Trans trans_cache_O24( legendre_cache_O24, Grid( "O24" ), 23, option::flt( false ) );
    }
}
#endif

CASE( "test_distspec" ) {
    trans::TransIFS trans( Grid( "F80" ), 159 );
    Log::info() << "Trans initialized" << std::endl;
    std::vector<double> rspecg;
    std::vector<int> nfrom;
    int nfld;
    Log::info() << "Read rspecg" << std::endl;
    read_rspecg( trans, rspecg, nfrom, nfld );

    std::vector<double> rspec( nfld * trans.trans()->nspec2 );
    std::vector<int> nto( nfld, 1 );
    std::vector<double> rgp( nfld * trans.trans()->ngptot );
    std::vector<double> rgpg( nfld * trans.grid().size() );
    std::vector<double> specnorms( nfld, 0 );

    trans.distspec( nfld, nfrom.data(), rspecg.data(), rspec.data() );
    trans.specnorm( nfld, rspec.data(), specnorms.data() );
    trans.invtrans( nfld, rspec.data(), rgp.data() );
    trans.gathgrid( nfld, nto.data(), rgp.data(), rgpg.data() );

    if ( mpi::comm().rank() == 0 ) {
        EXPECT( eckit::types::is_approximately_equal( specnorms[0], 1., 1.e-10 ) );
        EXPECT( eckit::types::is_approximately_equal( specnorms[1], 2., 1.e-10 ) );
    }

    Log::info() << "end test_distspec" << std::endl;
}

CASE( "test_distspec_speconly" ) {
    functionspace::Spectral fs( 159 );
    int nfld  = 2;
    Field glb = fs.createField<double>( option::global() | option::levels( nfld ) );
    read_rspecg( glb );

    std::vector<double> specnorms( nfld, 0 );

    Field dist = fs.createField( glb );

    fs.scatter( glb, dist );
    fs.norm( dist, specnorms );

    if ( mpi::comm().rank() == 0 ) {
        EXPECT( eckit::types::is_approximately_equal( specnorms[0], 1., 1.e-10 ) );
        EXPECT( eckit::types::is_approximately_equal( specnorms[1], 2., 1.e-10 ) );
    }
    Log::info() << "end test_distspec_only" << std::endl;
}

CASE( "test_distribution" ) {
    Grid g( "O80" );

    Log::info() << "test_distribution" << std::endl;

    grid::Distribution d_trans = grid::Partitioner( new TransPartitioner() ).partition( g );
    Log::info() << "trans distribution created" << std::endl;

    grid::Distribution d_eqreg = grid::Partitioner( new EqualRegionsPartitioner() ).partition( g );
    Log::info() << "eqregions distribution created" << std::endl;

    if ( mpi::comm().rank() == 0 ) {
        EXPECT( d_trans.nb_partitions() == d_eqreg.nb_partitions() );
        EXPECT( d_trans.max_pts() == d_eqreg.max_pts() );
        EXPECT( d_trans.min_pts() == d_eqreg.min_pts() );

        EXPECT( d_trans.nb_pts() == d_eqreg.nb_pts() );
    }
}

CASE( "test_generate_mesh" ) {
    Log::info() << "test_generate_mesh" << std::endl;
    Grid g( "O80" );
    StructuredMeshGenerator generate( atlas::util::Config( "angle", 0 )( "triangulate", true ) );

    Mesh m_default = generate( g );

    Log::info() << "trans_distribution" << std::endl;
    grid::Distribution trans_distribution = grid::Partitioner( new TransPartitioner() ).partition( g );
    Mesh m_trans                          = generate( g, trans_distribution );

    Log::info() << "eqreg_distribution" << std::endl;
    grid::Distribution eqreg_distribution = grid::Partitioner( new EqualRegionsPartitioner() ).partition( g );
    Mesh m_eqreg                          = generate( g, eqreg_distribution );

    array::ArrayView<int, 1> p_default = array::make_view<int, 1>( m_default.nodes().partition() );
    array::ArrayView<int, 1> p_trans   = array::make_view<int, 1>( m_trans.nodes().partition() );
    array::ArrayView<int, 1> p_eqreg   = array::make_view<int, 1>( m_eqreg.nodes().partition() );

    for ( idx_t j = 0; j < p_default.shape( 0 ); ++j ) {
        EXPECT( p_default( j ) == p_trans( j ) );
        EXPECT( p_default( j ) == p_eqreg( j ) );
    }

    output::Gmsh( "N16_trans.msh" ).write( m_trans );
}

CASE( "test_spectral_fields" ) {
    Log::info() << "test_spectral_fields" << std::endl;

    Grid g( "O48" );
    StructuredMeshGenerator generate( atlas::util::Config( "angle", 0 )( "triangulate", false ) );
    Mesh m = generate( g );

    trans::Trans trans( g, 47 );

    functionspace::NodeColumns nodal( m );
    functionspace::Spectral spectral( trans );

    Field spf = spectral.createField<double>( option::name( "spf" ) );
    Field gpf = nodal.createField<double>( option::name( "gpf" ) );

    array::make_view<double, 1>( gpf ).assign( 0 );

    EXPECT_NO_THROW( trans.dirtrans( gpf, spf ) );
    EXPECT_NO_THROW( trans.invtrans( spf, gpf ) );

    FieldSet gpfields;
    gpfields.add( gpf );
    FieldSet spfields;
    spfields.add( spf );

    EXPECT_NO_THROW( trans.dirtrans( gpfields, spfields ) );
    EXPECT_NO_THROW( trans.invtrans( spfields, gpfields ) );

    gpfields.add( gpf );
    EXPECT_THROWS_AS( trans.dirtrans( gpfields, spfields ), eckit::Exception );
}

CASE( "test_nomesh" ) {
    Log::info() << "test_spectral_fields" << std::endl;

    Grid g( "O48" );
    trans::Trans trans( g, 47 );

    functionspace::Spectral spectral( trans );
    functionspace::StructuredColumns gridpoints( g, grid::Partitioner( "trans" ) );

    Field spfg = spectral.createField<double>( option::name( "spf" ) | option::global() );
    Field spf  = spectral.createField<double>( option::name( "spf" ) );
    Field gpf  = gridpoints.createField<double>( option::name( "gpf" ) );
    Field gpfg = gridpoints.createField<double>( option::name( "gpf" ) | option::global() );

    array::ArrayView<double, 1> spg = array::make_view<double, 1>( spfg );

    spectral.parallel_for( option::global(), [&]( int real, int imag, int n, int m ) {
        spg( real ) = +m * spectral.truncation() + n;
        spg( imag ) = ( n == 0 ? 0 : -m * spectral.truncation() + n );
    } );

    EXPECT_NO_THROW( spectral.scatter( spfg, spf ) );

    array::ArrayView<double, 1> sp = array::make_view<double, 1>( spf );

    spectral.parallel_for( [&]( idx_t real, idx_t imag, int n, int m ) {
        EXPECT( int( sp( real ) ) == +m * spectral.truncation() + n );
        EXPECT( int( sp( imag ) ) == ( n == 0 ? 0 : -m * spectral.truncation() + n ) );

        sp( real ) = ( n == 0 ? 4. : 0. );
        sp( imag ) = 0.;
    } );

    EXPECT_NO_THROW( trans.invtrans( spf, gpf ) );

    EXPECT_NO_THROW( gridpoints.gather( gpf, gpfg ) );

    if ( mpi::comm().rank() == 0 ) {
        array::ArrayView<double, 1> gpg = array::make_view<double, 1>( gpfg );
        for ( idx_t jp = 0; jp < gpg.size(); ++jp ) {
            EXPECT( is_approximately_equal( gpg( jp ), 4., 0.001 ) );
            Log::debug() << "gpg(" << jp << ")   :   " << gpg( jp ) << std::endl;
        }
    }

    EXPECT_NO_THROW( gridpoints.scatter( gpfg, gpf ) );

    EXPECT_NO_THROW( trans.dirtrans( gpf, spf ) );

    EXPECT_NO_THROW( spectral.gather( spf, spfg ) );

    spectral.parallel_for( option::global(), [&]( idx_t real, idx_t imag, int n ) {
        EXPECT( is_approximately_equal( spg( real ), ( n == 0 ? 4. : 0. ), 0.001 ) );
        EXPECT( is_approximately_equal( spg( imag ), 0., 0.001 ) );
    } );
}

CASE( "test_trans_factory" ) {
    Log::info() << "test_trans_factory" << std::endl;

    trans::TransFactory::list( Log::info() );
    Log::info() << std::endl;

    functionspace::StructuredColumns gp( Grid( "O48" ), grid::Partitioner( "trans" ) );
    functionspace::Spectral sp( 47 );

    trans::Trans trans1 = trans::Trans( gp, sp );
    EXPECT( bool( trans1 ) == true );

    trans::Trans trans2 = trans::Trans( Grid( "O48" ), 47 );
    EXPECT( bool( trans2 ) == true );
}

CASE( "test_trans_using_grid" ) {
    Log::info() << "test_trans_using_grid" << std::endl;

    trans::Trans trans( Grid( "O48" ), 47 );

    functionspace::StructuredColumns gp( trans.grid(), grid::Partitioner( "trans" ) );
    functionspace::Spectral sp( trans.truncation() );

    Field spf = sp.createField<double>( option::name( "spf" ) );
    Field gpf = gp.createField<double>( option::name( "gpf" ) );

    array::make_view<double, 1>( gpf ).assign( 0 );

    EXPECT_NO_THROW( trans.dirtrans( gpf, spf ) );
    EXPECT_NO_THROW( trans.invtrans( spf, gpf ) );

    FieldSet gpfields;
    gpfields.add( gpf );
    FieldSet spfields;
    spfields.add( spf );

    EXPECT_NO_THROW( trans.dirtrans( gpfields, spfields ) );
    EXPECT_NO_THROW( trans.invtrans( spfields, gpfields ) );

    gpfields.add( gpf );
    EXPECT_THROWS_AS( trans.dirtrans( gpfields, spfields ), eckit::Exception );
}

CASE( "test_trans_using_functionspace_NodeColumns" ) {
    Log::info() << "test_trans_using_functionspace_NodeColumns" << std::endl;

    functionspace::NodeColumns gp( MeshGenerator( "structured" ).generate( Grid( "O48" ) ) );
    functionspace::Spectral sp( 47 );

    trans::Trans trans( gp, sp );

    Field spf = sp.createField<double>( option::name( "spf" ) );
    Field gpf = gp.createField<double>( option::name( "gpf" ) );

    array::make_view<double, 1>( gpf ).assign( 0 );

    EXPECT_NO_THROW( trans.dirtrans( gpf, spf ) );
    EXPECT_NO_THROW( trans.invtrans( spf, gpf ) );

    FieldSet gpfields;
    gpfields.add( gpf );
    FieldSet spfields;
    spfields.add( spf );

    EXPECT_NO_THROW( trans.dirtrans( gpfields, spfields ) );
    EXPECT_NO_THROW( trans.invtrans( spfields, gpfields ) );

    gpfields.add( gpf );
    EXPECT_THROWS_AS( trans.dirtrans( gpfields, spfields ), eckit::Exception );
}

CASE( "test_trans_using_functionspace_StructuredColumns" ) {
    Log::info() << "test_trans_using_functionspace_StructuredColumns" << std::endl;

    functionspace::StructuredColumns gp( Grid( "O48" ), grid::Partitioner( "trans" ) );
    functionspace::Spectral sp( 47 );

    trans::Trans trans( gp, sp );

    Field spf = sp.createField<double>( option::name( "spf" ) );
    Field gpf = gp.createField<double>( option::name( "gpf" ) );

    array::make_view<double, 1>( gpf ).assign( 0 );

    EXPECT_NO_THROW( trans.dirtrans( gpf, spf ) );
    EXPECT_NO_THROW( trans.invtrans( spf, gpf ) );

    FieldSet gpfields;
    gpfields.add( gpf );
    FieldSet spfields;
    spfields.add( spf );

    EXPECT_NO_THROW( trans.dirtrans( gpfields, spfields ) );
    EXPECT_NO_THROW( trans.invtrans( spfields, gpfields ) );

    gpfields.add( gpf );
    EXPECT_THROWS_AS( trans.dirtrans( gpfields, spfields ), eckit::Exception );
}

CASE( "test_trans_MIR_lonlat" ) {
    Log::info() << "test_trans_MIR_lonlat" << std::endl;

    Grid grid( "L48" );
    trans::Trans trans( grid, 47 );

    // global fields
    std::vector<double> spf( trans.spectralCoefficients(), 0. );
    std::vector<double> gpf( grid.size() );

    if ( mpi::comm().size() == 1 ) {
        EXPECT_NO_THROW( trans.invtrans( 1, spf.data(), gpf.data(), option::global() ) );

        EXPECT_NO_THROW( trans.dirtrans( 1, gpf.data(), spf.data(), option::global() ) );
    }
}

CASE( "test_trans_VorDivToUV" ) {
    int nfld = 1;                          // TODO: test for nfld>1
    std::vector<int> truncation_array{1};  // truncation_array{159,160,1279};
    for ( size_t i = 0; i < truncation_array.size(); ++i ) {
        int truncation = truncation_array[i];
        int nspec2     = ( truncation + 1 ) * ( truncation + 2 );

        Log::info() << "truncation = " << truncation << std::endl;

        std::vector<double> field_vor( nfld * nspec2, 0. );
        std::vector<double> field_div( nfld * nspec2, 0. );

        // TODO: initialise field_vor and field_div with something meaningful
        field_vor[2 * nfld] = 1.;
        Log::info() << "vor: " << std::endl;
        for ( int j = 0; j < nfld * nspec2; j++ )
            Log::info() << field_vor[j] << " ";
        Log::info() << std::endl;

        // With IFS
        if ( trans::VorDivToUVFactory::has( "ifs" ) ) {
            trans::VorDivToUV vordiv_to_UV( truncation, option::type( "ifs" ) );
            EXPECT( vordiv_to_UV.truncation() == truncation );

            std::vector<double> field_U( nfld * nspec2 );
            std::vector<double> field_V( nfld * nspec2 );

            vordiv_to_UV.execute( nspec2, nfld, field_vor.data(), field_div.data(), field_U.data(), field_V.data() );

            // TODO: do some meaningful checks
            Log::info() << "Trans library" << std::endl;
            Log::info() << "U: " << std::endl;
            for ( int j = 0; j < nfld * nspec2; j++ )
                Log::info() << field_U[j] << " ";
            Log::info() << std::endl;
        }

        // With Local
        {
            trans::VorDivToUV vordiv_to_UV( truncation, option::type( "local" ) );
            EXPECT( vordiv_to_UV.truncation() == truncation );

            std::vector<double> field_U( nfld * nspec2 );
            std::vector<double> field_V( nfld * nspec2 );

            vordiv_to_UV.execute( nspec2, nfld, field_vor.data(), field_div.data(), field_U.data(), field_V.data() );

            // TODO: do some meaningful checks
            Log::info() << "Local transform" << std::endl;
            Log::info() << "U: " << std::endl;
            for ( int j = 0; j < nfld * nspec2; j++ )
                Log::info() << field_U[j] << " ";
            Log::info() << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run<atlas::test::AtlasTransEnvironment>( argc, argv );
}
