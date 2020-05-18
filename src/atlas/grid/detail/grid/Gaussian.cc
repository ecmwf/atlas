/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Gaussian.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>

#include "eckit/utils/Translator.h"

#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/pl/classic_gaussian/PointsPerLatitude.h"

namespace atlas {
namespace grid {
namespace {  // anonymous

static eckit::Translator<std::string, int> to_int;

static Projection projection( const Grid::Config& grid ) {
    Grid::Config config;
    if ( grid.get( "projection", config ) ) {
        return Projection( config );
    }
    return Projection();
}

static Domain domain( const Grid::Config& grid ) {
    Grid::Config config;
    if ( grid.get( "domain", config ) ) {
        return Domain( config );
    }
    return Domain();
}

static Spacing yspace( const Grid::Config& grid ) {
    long N;
    grid.get( "N", N );

    Grid::Config config;
    config.set( "type", "gaussian" );
    config.set( "start", 90.0 );
    config.set( "end", -90.0 );
    config.set( "N", 2 * N );
    return Spacing( config );
}

template <typename vector_t>
static StructuredGrid::XSpace xspace( const vector_t& nx ) {
    return StructuredGrid::XSpace( {0., 360.}, nx, false );
}

//---------------------------------------------------------------------------------------------------------------------

static class classic_gaussian : public GridBuilder {
public:
    classic_gaussian() : GridBuilder( "classic_gaussian", {"^[Nn]([0-9]+)$"}, {"N<gauss>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "N<gauss>"
           << "Classic Gaussian grid";
    }

    const Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        int id;
        std::vector<std::string> matches;
        if ( match( name, matches, id ) ) {
            util::Config grid( config );
            int N = to_int( matches[0] );
            grid.set( "type", type() );
            grid.set( "N", N );
            return create( grid );
        }
        return nullptr;
    }

    const Grid::Implementation* create( const Grid::Config& config ) const override {
        long N;
        config.get( "N", N );
        std::vector<idx_t> nx( 2 * N );
        detail::pl::classic_gaussian::points_per_latitude_npole_spole( N, nx.data() );
        return new StructuredGrid::grid_t( "N" + std::to_string( N ), xspace( nx ), yspace( config ),
                                           projection( config ), domain( config ) );
    }

    void force_link() {}

} classic_gaussian_;

//---------------------------------------------------------------------------------------------------------------------

static class octahedral_gaussian : GridBuilder {
public:
    octahedral_gaussian() : GridBuilder( "octahedral_gaussian", {"^[Oo]([0-9]+)$"}, {"O<gauss>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "O<gauss>"
           << "Octahedral Gaussian grid";
    }

    const Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        int id;
        std::vector<std::string> matches;
        if ( match( name, matches, id ) ) {
            util::Config grid( config );
            int N = to_int( matches[0] );
            grid.set( "type", type() );
            grid.set( "N", N );
            return create( grid );
        }
        return nullptr;
    }

    const Grid::Implementation* create( const Grid::Config& config ) const override {
        long N;
        config.get( "N", N );

        long start = 20;
        config.get( "nx[0]", start );

        std::vector<idx_t> nx( 2 * N );
        for ( long j = 0; j < N; ++j ) {
            nx[j]             = start + 4 * j;
            nx[2 * N - 1 - j] = nx[j];
        }
        return new StructuredGrid::grid_t( "O" + std::to_string( N ), xspace( nx ), yspace( config ),
                                           projection( config ), domain( config ) );
    }

    void force_link() {}

} octahedral_gaussian_;

//---------------------------------------------------------------------------------------------------------------------

static class regular_gaussian : GridBuilder {
public:
    regular_gaussian() : GridBuilder( "regular_gaussian", {"^[Ff]([0-9]+)$"}, {"F<gauss>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "F<gauss>"
           << "Regular Gaussian grid";
    }

    const Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        int id;
        std::vector<std::string> matches;
        if ( match( name, matches, id ) ) {
            util::Config grid( config );
            int N = to_int( matches[0] );
            grid.set( "type", type() );
            grid.set( "N", N );
            return create( grid );
        }
        return nullptr;
    }

    const Grid::Implementation* create( const Grid::Config& config ) const override {
        long N;
        config.get( "N", N );
        std::vector<idx_t> nx( 2 * N, 4 * N );
        return new StructuredGrid::grid_t( "F" + std::to_string( N ), xspace( nx ), yspace( config ),
                                           projection( config ), domain( config ) );
    }

    void force_link() {}

} regular_gaussian_;

//---------------------------------------------------------------------------------------------------------------------

}  // anonymous namespace

namespace detail {
namespace grid {

void force_link_Gaussian() {
    regular_gaussian_.force_link();
    classic_gaussian_.force_link();
    octahedral_gaussian_.force_link();
}

StructuredGrid::grid_t* reduced_gaussian( const std::vector<long>& nx ) {
    Grid::Config yspace;
    yspace.set( "type", "gaussian" );
    yspace.set( "start", 90.0 );
    yspace.set( "end", -90.0 );
    yspace.set( "N", nx.size() );

    return new StructuredGrid::grid_t( xspace( nx ), Spacing( yspace ), Projection(), Domain() );
}

StructuredGrid::grid_t* reduced_gaussian( const std::vector<long>& nx, const Domain& domain ) {
    Grid::Config yspace;
    yspace.set( "type", "gaussian" );
    yspace.set( "start", 90.0 );
    yspace.set( "end", -90.0 );
    yspace.set( "N", nx.size() );

    return new StructuredGrid::grid_t( xspace( nx ), Spacing( yspace ), Projection(), domain );
}

StructuredGrid::grid_t* reduced_gaussian( const std::vector<long>& nx, const Projection& projection ) {
    Grid::Config yspace;
    yspace.set( "type", "gaussian" );
    yspace.set( "start", 90.0 );
    yspace.set( "end", -90.0 );
    yspace.set( "N", nx.size() );

    return new StructuredGrid::grid_t( xspace( nx ), Spacing( yspace ), projection, Domain() );
}

StructuredGrid::grid_t* reduced_gaussian( const std::vector<int>& nx ) {
    Grid::Config yspace;
    yspace.set( "type", "gaussian" );
    yspace.set( "start", 90.0 );
    yspace.set( "end", -90.0 );
    yspace.set( "N", nx.size() );

    return new StructuredGrid::grid_t( xspace( nx ), Spacing( yspace ), Projection(), Domain() );
}

StructuredGrid::grid_t* reduced_gaussian( const std::vector<int>& nx, const Domain& domain ) {
    Grid::Config yspace;
    yspace.set( "type", "gaussian" );
    yspace.set( "start", 90.0 );
    yspace.set( "end", -90.0 );
    yspace.set( "N", nx.size() );

    return new StructuredGrid::grid_t( xspace( nx ), Spacing( yspace ), Projection(), domain );
}

StructuredGrid::grid_t* reduced_gaussian( const std::vector<int>& nx, const Projection& projection ) {
    Grid::Config yspace;
    yspace.set( "type", "gaussian" );
    yspace.set( "start", 90.0 );
    yspace.set( "end", -90.0 );
    yspace.set( "N", nx.size() );

    return new StructuredGrid::grid_t( xspace( nx ), Spacing( yspace ), projection, Domain() );
}


template <typename Int>
inline std::vector<idx_t> idx_vector( Int nx, idx_t ny ) {
    std::vector<idx_t> _nx( ny );
    for ( idx_t j = 0; j < ny; ++j ) {
        _nx[j] = nx[j];
    }
    return _nx;
}

extern "C" {

StructuredGrid::grid_t* atlas__grid__reduced__ReducedGaussian_int( int nx[], long ny ) {
    return reduced_gaussian( idx_vector( nx, ny ) );
}

StructuredGrid::grid_t* atlas__grid__reduced__ReducedGaussian_long( long nx[], long ny ) {
    return reduced_gaussian( idx_vector( nx, ny ) );
}

StructuredGrid::grid_t* atlas__grid__reduced__ReducedGaussian_int_projection(
    int nx[], long ny, const Projection::Implementation* projection ) {
    return reduced_gaussian( idx_vector( nx, ny ), Projection( projection ) );
}

StructuredGrid::grid_t* atlas__grid__reduced__ReducedGaussian_long_projection(
    long nx[], long ny, const Projection::Implementation* projection ) {
    return reduced_gaussian( idx_vector( nx, ny ), Projection( projection ) );
}
}

}  // namespace grid
}  // namespace detail

}  // namespace grid
}  // namespace atlas
