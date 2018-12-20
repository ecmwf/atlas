/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace/StructuredColumns.h"

#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>

#include "eckit/utils/MD5.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/Statistics.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Checksum.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/detail/Cache.h"

#define IDX( i, j ) "(" << i << "," << j << ")"

#define REMOTE_IDX_BASE 0

namespace atlas {
namespace functionspace {
namespace detail {

namespace {

template <typename T>
array::LocalView<T, 3> make_leveled_view( const Field& field ) {
    using namespace array;
    if ( field.levels() ) {
        if ( field.variables() ) { return make_view<T, 3>( field ).slice( Range::all(), Range::all(), Range::all() ); }
        else {
            return make_view<T, 2>( field ).slice( Range::all(), Range::all(), Range::dummy() );
        }
    }
    else {
        if ( field.variables() ) {
            return make_view<T, 2>( field ).slice( Range::all(), Range::dummy(), Range::all() );
        }
        else {
            return make_view<T, 1>( field ).slice( Range::all(), Range::dummy(), Range::dummy() );
        }
    }
}

template <typename T>
std::string checksum_3d_field( const parallel::Checksum& checksum, const Field& field ) {
    array::LocalView<T, 3> values = make_leveled_view<T>( field );
    array::ArrayT<T> surface_field( values.shape( 0 ), values.shape( 2 ) );
    array::ArrayView<T, 2> surface = array::make_view<T, 2>( surface_field );
    const idx_t npts               = values.shape( 0 );
    atlas_omp_for( idx_t n = 0; n < npts; ++n ) {
        for ( idx_t j = 0; j < surface.shape( 1 ); ++j ) {
            surface( n, j ) = 0.;
            for ( idx_t l = 0; l < values.shape( 1 ); ++l )
                surface( n, j ) += values( n, l, j );
        }
    }
    return checksum.execute( surface.data(), surface_field.stride( 0 ) );
}

struct GridPoint {
public:
    idx_t i, j;
    idx_t r;

    GridPoint( idx_t _i, idx_t _j ) : i( _i ), j( _j ) {}
    GridPoint( idx_t _i, idx_t _j, idx_t _r ) : i( _i ), j( _j ), r( _r ) {}

    /* No longer used:

    bool operator<( const GridPoint& other ) const {
        if ( j < other.j ) return true;
        if ( j == other.j ) return i < other.i;
        return false;
    }

    bool operator==( const GridPoint& other ) const { return ( j == other.j && i == other.i ); }

    */
};


struct GridPointSet {
public:
    GridPointSet() = default;
    GridPointSet( size_t size ) { set.reserve( size ); }

    std::vector<GridPoint> set;
    bool insert( idx_t i, idx_t j ) {
        idx_t r = set.size();
        set.emplace_back( i, j, r );
        return true;
    }

    idx_t size() const { return static_cast<idx_t>( set.size() ); }

    using const_iterator = decltype( set )::const_iterator;

    const_iterator begin() const { return set.begin(); }
    const_iterator end() const { return set.end(); }
};

}  // namespace

class StructuredColumnsHaloExchangeCache : public util::Cache<std::string, parallel::HaloExchange>
//        , public mesh::detail::MeshObserver
{
private:
    using Base = util::Cache<std::string, parallel::HaloExchange>;
    StructuredColumnsHaloExchangeCache() : Base( "StructuredColumnsHaloExchangeCache" ) {}

public:
    static StructuredColumnsHaloExchangeCache& instance() {
        static StructuredColumnsHaloExchangeCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create( const detail::StructuredColumns& grid ) {
        creator_type creator = std::bind( &StructuredColumnsHaloExchangeCache::create, &grid );
        return Base::get_or_create( key( grid ), creator );
    }
    virtual void onGridDestruction( detail::StructuredColumns& grid ) { remove( key( grid ) ); }

private:
    static Base::key_type key( const detail::StructuredColumns& grid ) {
        std::ostringstream key;
        key << "grid[address=" << &grid << "]";
        return key.str();
    }

    static value_type* create( const detail::StructuredColumns* grid ) {
        //mesh.get()->attachObserver( instance() );

        value_type* value = new value_type();

        value->setup( array::make_view<int, 1>( grid->partition() ).data(),
                      array::make_view<idx_t, 1>( grid->remote_index() ).data(), REMOTE_IDX_BASE, grid->sizeHalo() );
        return value;
    }
    virtual ~StructuredColumnsHaloExchangeCache() {}
};

class StructuredColumnsGatherScatterCache : public util::Cache<std::string, parallel::GatherScatter>
//        , public mesh::detail::MeshObserver
{
private:
    using Base = util::Cache<std::string, parallel::GatherScatter>;
    StructuredColumnsGatherScatterCache() : Base( "StructuredColumnsGatherScatterCache" ) {}

public:
    static StructuredColumnsGatherScatterCache& instance() {
        static StructuredColumnsGatherScatterCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create( const detail::StructuredColumns& grid ) {
        creator_type creator = std::bind( &StructuredColumnsGatherScatterCache::create, &grid );
        return Base::get_or_create( key( grid ), creator );
    }
    virtual void onGridDestruction( detail::StructuredColumns& grid ) { remove( key( grid ) ); }

private:
    static Base::key_type key( const detail::StructuredColumns& grid ) {
        std::ostringstream key;
        key << "grid[address=" << &grid << "]";
        return key.str();
    }

    static value_type* create( const detail::StructuredColumns* grid ) {
        //mesh.get()->attachObserver( instance() );

        value_type* value = new value_type();

        value->setup( array::make_view<int, 1>( grid->partition() ).data(),
                      array::make_view<idx_t, 1>( grid->remote_index() ).data(), REMOTE_IDX_BASE,
                      array::make_view<gidx_t, 1>( grid->global_index() ).data(), grid->sizeOwned() );
        return value;
    }
    virtual ~StructuredColumnsGatherScatterCache() {}
};

class StructuredColumnsChecksumCache : public util::Cache<std::string, parallel::Checksum>
//        , public mesh::detail::MeshObserver
{
private:
    using Base = util::Cache<std::string, parallel::Checksum>;
    StructuredColumnsChecksumCache() : Base( "StructuredColumnsChecksumCache" ) {}

public:
    static StructuredColumnsChecksumCache& instance() {
        static StructuredColumnsChecksumCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create( const detail::StructuredColumns& grid ) {
        creator_type creator = std::bind( &StructuredColumnsChecksumCache::create, &grid );
        return Base::get_or_create( key( grid ), creator );
    }
    //virtual void onMeshDestruction( mesh::detail::MeshImpl& mesh ) { remove( key( mesh ) ); }

private:
    static Base::key_type key( const detail::StructuredColumns& grid ) {
        std::ostringstream key;
        key << "mesh[address=" << &grid << "]";
        return key.str();
    }

    static value_type* create( const detail::StructuredColumns* grid ) {
        //mesh.get()->attachObserver( instance() );
        value_type* value = new value_type();
        util::ObjectHandle<parallel::GatherScatter> gather(
            StructuredColumnsGatherScatterCache::instance().get_or_create( *grid ) );
        value->setup( gather );
        return value;
    }
    virtual ~StructuredColumnsChecksumCache() {}
};


const parallel::GatherScatter& StructuredColumns::gather() const {
    if ( gather_scatter_ ) return *gather_scatter_;
    gather_scatter_ = StructuredColumnsGatherScatterCache::instance().get_or_create( *this );
    return *gather_scatter_;
}

const parallel::GatherScatter& StructuredColumns::scatter() const {
    if ( gather_scatter_ ) return *gather_scatter_;
    gather_scatter_ = StructuredColumnsGatherScatterCache::instance().get_or_create( *this );
    return *gather_scatter_;
}

const parallel::Checksum& StructuredColumns::checksum() const {
    if ( checksum_ ) return *checksum_;
    checksum_ = StructuredColumnsChecksumCache::instance().get_or_create( *this );
    return *checksum_;
}

const parallel::HaloExchange& StructuredColumns::halo_exchange() const {
    if ( halo_exchange_ ) return *halo_exchange_;
    halo_exchange_ = StructuredColumnsHaloExchangeCache::instance().get_or_create( *this );
    return *halo_exchange_;
}

void StructuredColumns::set_field_metadata( const eckit::Configuration& config, Field& field ) const {
    field.set_functionspace( this );

    bool global( false );
    if ( config.get( "global", global ) ) {
        if ( global ) {
            idx_t owner( 0 );
            config.get( "owner", owner );
            field.metadata().set( "owner", owner );
        }
    }
    field.metadata().set( "global", global );

    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    field.set_levels( levels );

    idx_t variables( 0 );
    config.get( "variables", variables );
    field.set_variables( variables );
}

array::DataType StructuredColumns::config_datatype( const eckit::Configuration& config ) const {
    array::DataType::kind_t kind;
    if ( !config.get( "datatype", kind ) ) throw eckit::AssertionFailed( "datatype missing", Here() );
    return array::DataType( kind );
}

std::string StructuredColumns::config_name( const eckit::Configuration& config ) const {
    std::string name;
    config.get( "name", name );
    return name;
}

idx_t StructuredColumns::config_levels( const eckit::Configuration& config ) const {
    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    return levels;
}

array::ArrayShape StructuredColumns::config_shape( const eckit::Configuration& config ) const {
    array::ArrayShape shape;

    shape.push_back( config_size( config ) );

    idx_t levels( nb_levels_ );
    config.get( "levels", levels );
    if ( levels > 0 ) shape.push_back( levels );

    idx_t variables( 0 );
    config.get( "variables", variables );
    if ( variables > 0 ) shape.push_back( variables );

    return shape;
}
size_t StructuredColumns::Map2to1::footprint() const {
    size_t size = sizeof( *this );
    size += data_.size() * sizeof( decltype( data_ )::value_type );
    return size;
}

void StructuredColumns::Map2to1::print( std::ostream& out ) const {
    for ( idx_t j = j_min_; j <= j_max_; ++j ) {
        out << std::setw( 4 ) << j << " : ";
        for ( idx_t i = i_min_; i <= i_max_; ++i ) {
            idx_t v = operator()( i, j );
            if ( v == missing() )
                out << std::setw( 6 ) << "X";
            else
                out << std::setw( 6 ) << v;
        }
        out << '\n';
    }
}

void StructuredColumns::IndexRange::print( std::ostream& out ) const {
    for ( idx_t i = min_; i <= max_; ++i ) {
        idx_t v = operator()( i );
        if ( v == missing() )
            out << std::setw( 4 ) << "X";
        else
            out << std::setw( 4 ) << v;
    }
    out << '\n';
}

idx_t StructuredColumns::config_size( const eckit::Configuration& config ) const {
    idx_t size = size_halo_;
    bool global( false );
    if ( config.get( "global", global ) ) {
        if ( global ) {
            idx_t owner( 0 );
            config.get( "owner", owner );
            size = ( static_cast<idx_t>( mpi::comm().rank() ) == owner ? grid_->size() : 0 );
        }
    }
    return size;
}

std::string StructuredColumns::distribution() const {
    return distribution_;
}

void StructuredColumns::throw_outofbounds( idx_t i, idx_t j ) const {
    std::stringstream ss;
    if ( j < j_begin_halo() || j >= j_end_halo() ) {
        ss << "OutofBounds: j out of range! : (i,j) = (" << i << "," << j << ") --- Expected: " << j_begin_halo()
           << " <= j < " << j_end_halo();
    }
    if ( i < i_begin_halo( j ) || i >= i_end_halo( j ) ) {
        ss << "OutofBounds: i out of range! : (i,j) = (" << i << "," << j << ") --- Expected: " << i_begin_halo( j )
           << " <= i < " << i_end_halo( j );
    }
    throw eckit::Exception( ss.str(), Here() );
}

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
StructuredColumns::StructuredColumns( const Grid& grid, const eckit::Configuration& config ) :
    StructuredColumns::StructuredColumns( grid, grid::Partitioner(), config ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Partitioner& p,
                                      const eckit::Configuration& config ) :
    StructuredColumns( grid, Vertical( config ), p, config ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Distribution& distribution,
                                      const eckit::Configuration& config ) :
    StructuredColumns( grid, distribution, Vertical( config ), config ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Distribution& distribution,
                                      const Vertical& vertical, const eckit::Configuration& config ) :
    vertical_( vertical ),
    nb_levels_( vertical_.size() ),
    grid_( new grid::StructuredGrid( grid ) ) {
    setup( distribution, config );
}

StructuredColumns::StructuredColumns( const Grid& grid, const Vertical& vertical, const eckit::Configuration& config ) :
    StructuredColumns( grid, vertical, grid::Partitioner(), config ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const Vertical& vertical, const grid::Partitioner& p,
                                      const eckit::Configuration& config ) :
    vertical_( vertical ),
    nb_levels_( vertical_.size() ),
    grid_( new grid::StructuredGrid( grid ) ) {
    grid::Partitioner partitioner( p );
    if ( not partitioner ) {
        if ( grid_->domain().global() ) { partitioner = grid::Partitioner( "equal_regions" ); }
        else {
            partitioner = grid::Partitioner( "checkerboard" );
        }
    }

    grid::Distribution distribution;
    ATLAS_TRACE_SCOPE( "Partitioning grid ..." ) { distribution = grid::Distribution( grid, partitioner ); }

    setup( distribution, config );
}

void StructuredColumns::setup( const grid::Distribution& distribution, const eckit::Configuration& config ) {
    ATLAS_TRACE( "Generating StructuredColumns..." );
    bool periodic_points = config.getInt( "periodic_points", false );
    if ( not( *grid_ ) ) { throw eckit::BadCast( "Grid is not a grid::Structured type", Here() ); }
    const eckit::mpi::Comm& comm = mpi::comm();


    ny_                  = grid_->ny();
    north_pole_included_ = 90. - grid_->y( 0 ) == 0.;
    south_pole_included_ = 90. + grid_->y( ny_ - 1 ) == 0;

    distribution_ = distribution.type();

    int mpi_rank = int( mpi::comm().rank() );

    j_begin_ = std::numeric_limits<idx_t>::max();
    j_end_   = std::numeric_limits<idx_t>::min();
    i_begin_.resize( grid_->ny(), std::numeric_limits<idx_t>::max() );
    i_end_.resize( grid_->ny(), std::numeric_limits<idx_t>::min() );
    idx_t c( 0 );
    idx_t owned( 0 );
    for ( idx_t j = 0; j < grid_->ny(); ++j ) {
        for ( idx_t i = 0; i < grid_->nx( j ); ++i, ++c ) {
            if ( distribution.partition( c ) == mpi_rank ) {
                j_begin_    = std::min<idx_t>( j_begin_, j );
                j_end_      = std::max<idx_t>( j_end_, j + 1 );
                i_begin_[j] = std::min<idx_t>( i_begin_[j], i );
                i_end_[j]   = std::max<idx_t>( i_end_[j], i + 1 );
                ++owned;
            }
        }
    }

    size_owned_ = owned;

    int halo = config.getInt( "halo", 0 );
    halo_    = halo;

    j_begin_halo_ = j_begin_ - halo;
    j_end_halo_   = j_end_ + halo;
    i_begin_halo_.resize( -halo, grid_->ny() - 1 + halo );
    i_end_halo_.resize( -halo, grid_->ny() - 1 + halo );

    auto compute_i = [this]( idx_t i, idx_t j ) -> idx_t {
        const idx_t nx = grid_->nx( j );
        while ( i >= nx ) {
            i -= nx;
        }
        while ( i < 0 ) {
            i += nx;
        }
        return i;
    };

    std::function<idx_t( idx_t )> compute_j;
    compute_j = [this, &compute_j]( idx_t j ) -> idx_t {
        if ( j < 0 ) { j = ( grid_->y( 0 ) == 90. ) ? -j : -j - 1; }
        else if ( j >= grid_->ny() ) {
            idx_t jlast = grid_->ny() - 1;
            j           = ( grid_->y( jlast ) == -90. ) ? jlast - 1 - ( j - grid_->ny() ) : jlast - ( j - grid_->ny() );
        }
        if ( j < 0 or j >= grid_->ny() ) { j = compute_j( j ); }
        return j;
    };

    auto compute_x = [this, &compute_i, &compute_j]( idx_t i, idx_t j ) -> double {
        idx_t ii, jj;
        double x;
        jj             = compute_j( j );
        ii             = compute_i( i, jj );  // guaranteed between 0 and nx(jj)
        const idx_t nx = grid_->nx( jj );
        const double a = ( ii - i ) / nx;
        x              = grid_->x( ii, jj ) - a * grid_->x( nx, jj );
        return x;
    };

    auto compute_y = [this, &compute_j]( idx_t j ) -> double {
        idx_t jj;
        double y;
        jj = compute_j( j );
        y  = ( j < 0 ) ? 90. + ( 90. - grid_->y( jj ) )
                      : ( j >= grid_->ny() ) ? -90. + ( -90. - grid_->y( jj ) ) : grid_->y( jj );
        return y;
    };

    std::vector<gidx_t> global_offsets( grid_->ny() );
    idx_t grid_idx = 0;
    for ( idx_t j = 0; j < grid_->ny(); ++j ) {
        global_offsets[j] = grid_idx;
        grid_idx += grid_->nx( j );
    }

    auto compute_g = [this, &global_offsets, &compute_i, &compute_j]( idx_t i, idx_t j ) -> gidx_t {
        idx_t ii, jj;
        gidx_t g;
        jj = compute_j( j );
        ii = compute_i( i, jj );
        if ( jj != j ) {
            ASSERT( grid_->nx( jj ) % 2 == 0 );  // assert even number of points
            ii = ( ii < grid_->nx( jj ) / 2 ) ? ii + grid_->nx( jj ) / 2
                                              : ( ii >= grid_->nx( jj ) / 2 ) ? ii - grid_->nx( jj ) / 2 : ii;
        }
        g = global_offsets[jj] + ii + 1;
        return g;
    };

    auto compute_p = [this, &global_offsets, &distribution, &compute_i, &compute_j]( idx_t i, idx_t j ) -> int {
        idx_t ii, jj;
        int p;
        jj = compute_j( j );
        ii = compute_i( i, jj );
        if ( jj != j ) {
            ASSERT( grid_->nx( jj ) % 2 == 0 );  // assert even number of points
            ii = ( ii < grid_->nx( jj ) / 2 ) ? ii + grid_->nx( jj ) / 2
                                              : ( ii >= grid_->nx( jj ) / 2 ) ? ii - grid_->nx( jj ) / 2 : ii;
        }
        p = distribution.partition( global_offsets[jj] + ii );
        return p;
    };

    if ( atlas::Library::instance().debug() ) {
        ATLAS_TRACE_SCOPE( "Load imbalance" ) { comm.barrier(); }
    }

    GridPointSet gridpoints;

    ATLAS_TRACE_SCOPE( "Compute mapping ..." ) {
        idx_t imin = std::numeric_limits<idx_t>::max();
        idx_t imax = -std::numeric_limits<idx_t>::max();
        idx_t jmin = std::numeric_limits<idx_t>::max();
        idx_t jmax = -std::numeric_limits<idx_t>::max();

        ATLAS_TRACE_SCOPE( "Compute bounds" ) {
            for ( idx_t j = j_begin_halo_; j < j_end_halo_; ++j ) {
                i_begin_halo_( j ) = imin;
                i_end_halo_( j )   = imax;
            }
            double eps = 1.e-12;
            for ( idx_t j = j_begin_; j < j_end_; ++j ) {
                for ( idx_t i : {i_begin_[j], i_end_[j] - 1} ) {
                    // Following line only, increases periodic halo on the east side by 1
                    if ( periodic_points && i == grid_->nx( j ) - 1 ) { ++i; }

                    double x      = grid_->x( i, j );
                    double x_next = grid_->x( i + 1, j );
                    double x_prev = grid_->x( i - 1, j );
                    for ( idx_t jj = j - halo; jj <= j + halo; ++jj ) {
                        idx_t last = grid_->nx( compute_j( jj ) ) - 1;
                        if ( i == grid_->nx( j ) ) { ++last; }

                        jmin = std::min( jmin, jj );
                        jmax = std::max( jmax, jj );
                        // Compute ii as index less-equal of x
                        //
                        //              x(i,j)
                        //    |-----|-----|-----|-----|
                        // ii-halo       ii
                        //
                        //                 x(i,j)
                        //    |-----|-----|--+--|-----|
                        // ii-halo       ii

                        idx_t ii = -halo;
                        while ( compute_x( ii, jj ) < x - eps ) {
                            ii++;
                        }

                        // ATLAS-186 workaround
                        // This while should not have to be there, but is here because of
                        // the MatchingMeshDomainDecomposition algorithm. that may lead to points
                        // left of the point ii.
                        while ( compute_x( ii - 1, jj ) > x_prev + eps ) {
                            --ii;
                        }

                        idx_t i_minus_halo = ii - halo;

                        // Compute ii as index less-equal of x_next
                        //
                        //               x(i,j) x_next(i,j)
                        //   |-----|-----|-+---|-+---|-----|
                        // ii-halo            iii       iii+halo
                        //
                        idx_t iii = ii;
                        while ( compute_x( iii + 1, jj ) < x_next - eps ) {
                            ++iii;
                        }
                        iii               = std::min( iii, last );
                        idx_t i_plus_halo = iii + halo;

                        imin                = std::min( imin, i_minus_halo );
                        imax                = std::max( imax, i_plus_halo );
                        i_begin_halo_( jj ) = std::min( i_begin_halo_( jj ), i_minus_halo );
                        i_end_halo_( jj )   = std::max( i_end_halo_( jj ), i_plus_halo + 1 );
                    }
                }
            }
        }

        int extra_halo{0};
        for ( idx_t j = j_begin_halo_; j < j_begin_; ++j ) {
            extra_halo += i_end_halo_( j ) - i_begin_halo_( j );
        }
        for ( idx_t j = j_begin_; j < j_end_; ++j ) {
            extra_halo += i_begin_[j] - i_begin_halo_( j );
            extra_halo += i_end_halo_( j ) - i_end_[j];
        }
        for ( idx_t j = j_end_; j < j_end_halo_; ++j ) {
            extra_halo += i_end_halo_( j ) - i_begin_halo_( j );
        }


        gridpoints = GridPointSet{static_cast<size_t>( owned + extra_halo )};

        ATLAS_TRACE_SCOPE( "Assemble gridpoints" ) {
            for ( idx_t j = j_begin_; j < j_end_; ++j ) {
                for ( idx_t i = i_begin_[j]; i < i_end_[j]; ++i ) {
                    gridpoints.insert( i, j );
                }
            }

            ASSERT( gridpoints.size() == owned );

            for ( idx_t j = j_begin_halo_; j < j_begin_; ++j ) {
                for ( idx_t i = i_begin_halo_( j ); i < i_end_halo_( j ); ++i ) {
                    gridpoints.insert( i, j );
                }
            }
            for ( idx_t j = j_begin_; j < j_end_; ++j ) {
                for ( idx_t i = i_begin_halo_( j ); i < i_begin_[j]; ++i ) {
                    gridpoints.insert( i, j );
                }
                for ( idx_t i = i_end_[j]; i < i_end_halo_( j ); ++i ) {
                    gridpoints.insert( i, j );
                }
            }
            for ( idx_t j = j_end_; j < j_end_halo_; ++j ) {
                for ( idx_t i = i_begin_halo_( j ); i < i_end_halo_( j ); ++i ) {
                    gridpoints.insert( i, j );
                }
            }

            ASSERT( gridpoints.size() == owned + extra_halo );
        }

        ATLAS_TRACE_SCOPE( "Fill in ij2gp " ) {
            ij2gp_.resize( {imin, imax}, {jmin, jmax} );

            for ( const GridPoint& gp : gridpoints ) {
                ij2gp_.set( gp.i, gp.j, gp.r );
            }
        }
    }

    ATLAS_TRACE_SCOPE( "Create fields" ) {
        size_halo_          = gridpoints.size();
        field_partition_    = Field( "partition", array::make_datatype<int>(), array::make_shape( size_halo_ ) );
        field_global_index_ = Field( "glb_idx", array::make_datatype<gidx_t>(), array::make_shape( size_halo_ ) );
        field_index_i_      = Field( "index_i", array::make_datatype<idx_t>(), array::make_shape( size_halo_ ) );
        field_index_j_      = Field( "index_j", array::make_datatype<idx_t>(), array::make_shape( size_halo_ ) );
        field_xy_           = Field( "xy", array::make_datatype<double>(), array::make_shape( size_halo_, 2 ) );

        auto xy         = array::make_view<double, 2>( field_xy_ );
        auto part       = array::make_view<int, 1>( field_partition_ );
        auto global_idx = array::make_view<gidx_t, 1>( field_global_index_ );
        auto index_i    = array::make_indexview<idx_t, 1>( field_index_i_ );
        auto index_j    = array::make_indexview<idx_t, 1>( field_index_j_ );

        for ( const GridPoint& gp : gridpoints ) {
            xy( gp.r, XX ) = compute_x( gp.i, gp.j );
            if ( gp.j >= 0 && gp.j < grid_->ny() ) { xy( gp.r, YY ) = grid_->y( gp.j ); }
            else {
                xy( gp.r, YY ) = compute_y( gp.j );
            }

            bool in_domain( false );
            if ( gp.j >= 0 && gp.j < grid_->ny() ) {
                if ( gp.i >= 0 && gp.i < grid_->nx( gp.j ) ) {
                    in_domain          = true;
                    gidx_t k           = global_offsets[gp.j] + gp.i;
                    part( gp.r )       = distribution.partition( k );
                    global_idx( gp.r ) = k + 1;
                }
            }
            if ( not in_domain ) {
                global_idx( gp.r ) = compute_g( gp.i, gp.j );
                part( gp.r )       = compute_p( gp.i, gp.j );
            }
            index_i( gp.r ) = gp.i;
            index_j( gp.r ) = gp.j;
        }
    }
}


void StructuredColumns::create_remote_index() const {
    field_remote_index_ = Field( "remote_idx", array::make_datatype<idx_t>(), array::make_shape( size_halo_ ) );
    auto remote_idx     = array::make_view<idx_t, 1>( field_remote_index_ );
    for ( idx_t n = 0; n < size_owned_; ++n ) {
        remote_idx( n ) = n;
    }


    ATLAS_TRACE_SCOPE( "Parallelisation ..." ) {
        auto build_partition_graph = [this]() -> std::unique_ptr<Mesh::PartitionGraph> {
            const eckit::mpi::Comm& comm = mpi::comm();
            const int mpi_size           = int( comm.size() );
            const int mpi_rank           = int( comm.rank() );

            auto p = array::make_view<int, 1>( this->partition() );

            std::set<int> others_set;
            others_set.insert( mpi_rank );
            for ( idx_t i = size_owned_; i < size_halo_; ++i ) {
                others_set.insert( p( i ) );
            }
            std::vector<int> others( others_set.begin(), others_set.end() );

            eckit::mpi::Buffer<int> recv_others( mpi_size );

            comm.allGatherv( others.begin(), others.end(), recv_others );

            std::vector<idx_t> counts( recv_others.counts.begin(), recv_others.counts.end() );
            std::vector<idx_t> displs( recv_others.displs.begin(), recv_others.displs.end() );
            std::vector<idx_t> values( recv_others.buffer.begin(), recv_others.buffer.end() );
            return std::unique_ptr<Mesh::PartitionGraph>(
                new Mesh::PartitionGraph( values.data(), mpi_size, displs.data(), counts.data() ) );
        };

        std::unique_ptr<Mesh::PartitionGraph> graph_ptr;
        ATLAS_TRACE_SCOPE( "Building partition graph..." ) { graph_ptr = build_partition_graph(); }
        const Mesh::PartitionGraph& graph = *graph_ptr;

        ATLAS_TRACE_SCOPE( "Setup remote_index fields..." ) {
            auto p = array::make_view<int, 1>( partition() );
            auto g = array::make_view<gidx_t, 1>( global_index() );

            const eckit::mpi::Comm& comm = mpi::comm();
            const int mpi_rank           = int( comm.rank() );

            auto neighbours           = graph.nearestNeighbours( mpi_rank );
            const idx_t nb_neighbours = static_cast<idx_t>( neighbours.size() );
            std::unordered_map<int, idx_t> part_to_neighbour;
            part_to_neighbour.reserve( nb_neighbours );
            ATLAS_TRACE_SCOPE( "part_to_neighbour" )
            for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                part_to_neighbour[neighbours[j]] = j;
            }
            std::vector<idx_t> halo_per_neighbour( neighbours.size(), 0 );
            ATLAS_TRACE_SCOPE( "set halo_per_neighbour" )
            for ( idx_t i = size_owned_; i < size_halo_; ++i ) {
                halo_per_neighbour[part_to_neighbour[p( i )]]++;
            }

            std::vector<std::vector<gidx_t>> g_per_neighbour( neighbours.size() );
            for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                g_per_neighbour[j].reserve( halo_per_neighbour[j] );
            }
            for ( idx_t j = size_owned_; j < size_halo_; ++j ) {
                g_per_neighbour[part_to_neighbour[p( j )]].emplace_back( g( j ) );
            }
            std::vector<std::vector<idx_t>> r_per_neighbour( neighbours.size() );
            for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                r_per_neighbour[j].resize( halo_per_neighbour[j] );
            }

            std::vector<eckit::mpi::Request> send_requests( neighbours.size() );
            std::vector<eckit::mpi::Request> recv_requests( neighbours.size() );

            std::vector<idx_t> recv_size( neighbours.size() );
            int tag = 0;
            ATLAS_TRACE_SCOPE( "send-receive g_per_neighbour size" )
            for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                idx_t g_per_neighbour_size = static_cast<idx_t>( g_per_neighbour[j].size() );
                send_requests[j]           = comm.iSend( g_per_neighbour_size, neighbours[j], tag );
                recv_requests[j]           = comm.iReceive( recv_size[j], neighbours[j], tag );
            }

            ATLAS_TRACE_MPI( WAIT ) {
                for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                    comm.wait( send_requests[j] );
                }

                for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                    comm.wait( recv_requests[j] );
                }
            }

            std::vector<std::vector<gidx_t>> recv_g_per_neighbour( neighbours.size() );
            ATLAS_TRACE_SCOPE( "send-receive g_per_neighbour" )
            for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                recv_g_per_neighbour[j].resize( recv_size[j] );

                send_requests[j] =
                    comm.iSend( g_per_neighbour[j].data(), g_per_neighbour[j].size(), neighbours[j], tag );
                recv_requests[j] =
                    comm.iReceive( recv_g_per_neighbour[j].data(), recv_g_per_neighbour[j].size(), neighbours[j], tag );
            }

            std::vector<std::vector<idx_t>> send_r_per_neighbour( neighbours.size() );

            // Assemble "send_r_per_neighbour"
            // Depending if we can afford memory for a globally sized array, we can have a
            // much faster version of g_to_r map using std::vector.
            // TODO: using c++14 we can create a polymorphic lambda to avoid duplicated
            // code in the two branches of following if.
            if ( comm.size() == 1 ) {
                std::vector<idx_t> g_to_r( size_owned_ + 1 );

                ATLAS_TRACE_SCOPE( "g_to_r (using vector)" )
                for ( idx_t j = 0; j < size_owned_; ++j ) {
#if ATLAS_ARRAYVIEW_BOUNDS_CHECKING
                    if ( g( j ) >= size_owned_ + 1 ) {
                        ATLAS_DEBUG_VAR( g( j ) );
                        throw eckit::OutOfRange( g( j ), size_owned_ + 1, Here() );
                    }
#endif
                    g_to_r[g( j )] = j;
                }
                for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                    send_r_per_neighbour[j].reserve( recv_size[j] );

                    comm.wait( recv_requests[j] );  // wait for recv_g_per_neighbour[j]
                    for ( gidx_t gidx : recv_g_per_neighbour[j] ) {
                        send_r_per_neighbour[j].emplace_back( g_to_r[gidx] );
                    }
                }
            }
            else {
                std::unordered_map<gidx_t, idx_t> g_to_r;
                g_to_r.reserve( size_owned_ );

                ATLAS_TRACE_SCOPE( "g_to_r (using unordered set)" )
                for ( idx_t j = 0; j < size_owned_; ++j ) {
                    g_to_r[g( j )] = j;
                }
                for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                    send_r_per_neighbour[j].reserve( recv_size[j] );

                    comm.wait( recv_requests[j] );  // wait for recv_g_per_neighbour[j]
                    for ( gidx_t gidx : recv_g_per_neighbour[j] ) {
                        send_r_per_neighbour[j].emplace_back( g_to_r[gidx] );
                    }
                }
            }

            ATLAS_TRACE_MPI( ALLTOALL ) {
                for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                    comm.wait( send_requests[j] );
                    send_requests[j] = comm.iSend( send_r_per_neighbour[j].data(), send_r_per_neighbour[j].size(),
                                                   neighbours[j], tag );
                    recv_requests[j] =
                        comm.iReceive( r_per_neighbour[j].data(), r_per_neighbour[j].size(), neighbours[j], tag );
                }
            }

            ATLAS_TRACE_MPI( WAIT ) {
                for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                    comm.wait( recv_requests[j] );
                }
            }

            std::vector<idx_t> counters( neighbours.size(), 0 );
            ATLAS_TRACE_SCOPE( "set remote_idx" )
            for ( idx_t j = size_owned_; j < size_halo_; ++j ) {
                idx_t neighbour = part_to_neighbour[p( j )];
                remote_idx( j ) = r_per_neighbour[neighbour][counters[neighbour]++];
            }

            ATLAS_TRACE_MPI( WAIT ) {
                for ( idx_t j = 0; j < nb_neighbours; ++j ) {
                    comm.wait( send_requests[j] );
                }
            }
        }
    }
}


void StructuredColumns::compute_xy( idx_t i, idx_t j, PointXY& xy ) const {
    idx_t jj;
    if ( j < 0 ) {
        jj     = -j - 1 + north_pole_included_;
        xy.y() = 180. - grid_->y( jj );
    }
    else if ( j >= ny_ ) {
        jj     = 2 * ny_ - j - 1 - south_pole_included_;
        xy.y() = -180. - grid_->y( jj );
    }
    else {
        jj     = j;
        xy.y() = grid_->y( jj );
    }
    xy.x() = grid_->x( i, jj );
}


// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------
StructuredColumns::~StructuredColumns() {
    delete grid_;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Create Field
// ----------------------------------------------------------------------------
Field StructuredColumns::createField( const eckit::Configuration& options ) const {
    Field field( config_name( options ), config_datatype( options ), config_shape( options ) );
    set_field_metadata( options, field );
    return field;
}

Field StructuredColumns::createField( const Field& other, const eckit::Configuration& config ) const {
    return createField( option::datatype( other.datatype() ) | option::levels( other.levels() ) |
                        option::variables( other.variables() ) | config );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Gather FieldSet
// ----------------------------------------------------------------------------
void StructuredColumns::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const {
    ASSERT( local_fieldset.size() == global_fieldset.size() );

    for ( idx_t f = 0; f < local_fieldset.size(); ++f ) {
        const Field& loc      = local_fieldset[f];
        Field& glb            = global_fieldset[f];
        const idx_t nb_fields = 1;
        idx_t root( 0 );
        glb.metadata().get( "owner", root );

        if ( loc.datatype() == array::DataType::kind<int>() ) {
            parallel::Field<int const> loc_field( make_leveled_view<int>( loc ) );
            parallel::Field<int> glb_field( make_leveled_view<int>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<long>() ) {
            parallel::Field<long const> loc_field( make_leveled_view<long>( loc ) );
            parallel::Field<long> glb_field( make_leveled_view<long>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<float>() ) {
            parallel::Field<float const> loc_field( make_leveled_view<float>( loc ) );
            parallel::Field<float> glb_field( make_leveled_view<float>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<double>() ) {
            parallel::Field<double const> loc_field( make_leveled_view<double>( loc ) );
            parallel::Field<double> glb_field( make_leveled_view<double>( glb ) );
            gather().gather( &loc_field, &glb_field, nb_fields, root );
        }
        else
            throw eckit::Exception( "datatype not supported", Here() );
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Gather Field
// ----------------------------------------------------------------------------
void StructuredColumns::gather( const Field& local, Field& global ) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add( local );
    global_fields.add( global );
    gather( local_fields, global_fields );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Scatter FieldSet
// ----------------------------------------------------------------------------
void StructuredColumns::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const {
    ASSERT( local_fieldset.size() == global_fieldset.size() );

    for ( idx_t f = 0; f < local_fieldset.size(); ++f ) {
        const Field& glb      = global_fieldset[f];
        Field& loc            = local_fieldset[f];
        const idx_t nb_fields = 1;
        idx_t root( 0 );
        glb.metadata().get( "owner", root );

        if ( loc.datatype() == array::DataType::kind<int>() ) {
            parallel::Field<int const> glb_field( make_leveled_view<int>( glb ) );
            parallel::Field<int> loc_field( make_leveled_view<int>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<long>() ) {
            parallel::Field<long const> glb_field( make_leveled_view<long>( glb ) );
            parallel::Field<long> loc_field( make_leveled_view<long>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<float>() ) {
            parallel::Field<float const> glb_field( make_leveled_view<float>( glb ) );
            parallel::Field<float> loc_field( make_leveled_view<float>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else if ( loc.datatype() == array::DataType::kind<double>() ) {
            parallel::Field<double const> glb_field( make_leveled_view<double>( glb ) );
            parallel::Field<double> loc_field( make_leveled_view<double>( loc ) );
            scatter().scatter( &glb_field, &loc_field, nb_fields, root );
        }
        else
            throw eckit::Exception( "datatype not supported", Here() );

        glb.metadata().broadcast( loc.metadata(), root );
        loc.metadata().set( "global", false );
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Scatter Field
// ----------------------------------------------------------------------------
void StructuredColumns::scatter( const Field& global, Field& local ) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add( global );
    local_fields.add( local );
    scatter( global_fields, local_fields );
}
// ----------------------------------------------------------------------------

std::string StructuredColumns::checksum( const FieldSet& fieldset ) const {
    eckit::MD5 md5;
    for ( idx_t f = 0; f < fieldset.size(); ++f ) {
        const Field& field = fieldset[f];
        if ( field.datatype() == array::DataType::kind<int>() )
            md5 << checksum_3d_field<int>( checksum(), field );
        else if ( field.datatype() == array::DataType::kind<long>() )
            md5 << checksum_3d_field<long>( checksum(), field );
        else if ( field.datatype() == array::DataType::kind<float>() )
            md5 << checksum_3d_field<float>( checksum(), field );
        else if ( field.datatype() == array::DataType::kind<double>() )
            md5 << checksum_3d_field<double>( checksum(), field );
        else
            throw eckit::Exception( "datatype not supported", Here() );
    }
    return md5;
}
std::string StructuredColumns::checksum( const Field& field ) const {
    FieldSet fieldset;
    fieldset.add( field );
    return checksum( fieldset );
}

const grid::StructuredGrid& StructuredColumns::grid() const {
    return *grid_;
}

namespace {


template <int RANK>
struct FixupHaloForVectors {
    FixupHaloForVectors( const StructuredColumns& ) {}
    template <typename DATATYPE>
    void apply( Field& field ) {
        std::string type = field.metadata().getString( "type", "scalar" );
        if ( type == "vector " ) { NOTIMP; }
    }
};

template <>
struct FixupHaloForVectors<2> {
    static constexpr int RANK = 2;
    const StructuredColumns& fs;
    FixupHaloForVectors( const StructuredColumns& _fs ) : fs( _fs ) {}

    template <typename DATATYPE>
    void apply( Field& field ) {
        std::string type = field.metadata().getString( "type", "scalar" );
        if ( type == "vector" ) {
            auto array = array::make_view<DATATYPE, RANK>( field );
            for ( idx_t j = fs.j_begin_halo(); j < 0; ++j ) {
                for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                    idx_t n        = fs.index( i, j );
                    array( n, XX ) = -array( n, XX );
                    array( n, YY ) = -array( n, YY );
                }
            }
            for ( idx_t j = fs.grid().ny(); j < fs.j_end_halo(); ++j ) {
                for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                    idx_t n        = fs.index( i, j );
                    array( n, XX ) = -array( n, XX );
                    array( n, YY ) = -array( n, YY );
                }
            }
        }
    }
};

template <>
struct FixupHaloForVectors<3> {
    static constexpr int RANK = 3;
    const StructuredColumns& fs;
    FixupHaloForVectors( const StructuredColumns& _fs ) : fs( _fs ) {}

    template <typename DATATYPE>
    void apply( Field& field ) {
        std::string type = field.metadata().getString( "type", "scalar" );
        if ( type == "vector" ) {
            auto array = array::make_view<DATATYPE, RANK>( field );
            for ( idx_t j = fs.j_begin_halo(); j < 0; ++j ) {
                for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                    idx_t n = fs.index( i, j );
                    for ( idx_t k = fs.k_begin(); k < fs.k_end(); ++k ) {
                        array( n, k, XX ) = -array( n, k, XX );
                        array( n, k, YY ) = -array( n, k, YY );
                    }
                }
            }
            for ( idx_t j = fs.grid().ny(); j < fs.j_end_halo(); ++j ) {
                for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                    idx_t n = fs.index( i, j );
                    for ( idx_t k = fs.k_begin(); k < fs.k_end(); ++k ) {
                        array( n, k, XX ) = -array( n, k, XX );
                        array( n, k, YY ) = -array( n, k, YY );
                    }
                }
            }
        }
    }
};


template <int RANK>
void dispatch_haloExchange( Field& field, const parallel::HaloExchange& halo_exchange, const StructuredColumns& fs ) {
    FixupHaloForVectors<RANK> fixup_halos( fs );
    if ( field.datatype() == array::DataType::kind<int>() ) {
        halo_exchange.template execute<int, RANK>( field.array(), false );
        fixup_halos.template apply<int>( field );
    }
    else if ( field.datatype() == array::DataType::kind<long>() ) {
        halo_exchange.template execute<long, RANK>( field.array(), false );
        fixup_halos.template apply<long>( field );
    }
    else if ( field.datatype() == array::DataType::kind<float>() ) {
        halo_exchange.template execute<float, RANK>( field.array(), false );
        fixup_halos.template apply<float>( field );
    }
    else if ( field.datatype() == array::DataType::kind<double>() ) {
        halo_exchange.template execute<double, RANK>( field.array(), false );
        fixup_halos.template apply<double>( field );
    }
    else
        throw eckit::Exception( "datatype not supported", Here() );
    field.set_dirty( false );
}
}  // namespace

void StructuredColumns::haloExchange( const FieldSet& fieldset, bool ) const {
    for ( idx_t f = 0; f < fieldset.size(); ++f ) {
        Field& field = const_cast<FieldSet&>( fieldset )[f];
        switch ( field.rank() ) {
            case 1:
                dispatch_haloExchange<1>( field, halo_exchange(), *this );
                break;
            case 2:
                dispatch_haloExchange<2>( field, halo_exchange(), *this );
                break;
            case 3:
                dispatch_haloExchange<3>( field, halo_exchange(), *this );
                break;
            case 4:
                dispatch_haloExchange<4>( field, halo_exchange(), *this );
                break;
            default:
                throw eckit::Exception( "Rank not supported", Here() );
        }
    }
}

void StructuredColumns::haloExchange( const Field& field, bool ) const {
    FieldSet fieldset;
    fieldset.add( field );
    haloExchange( fieldset );
}

size_t StructuredColumns::footprint() const {
    size_t size = sizeof( *this );
    size += ij2gp_.footprint();
    if ( field_xy_ ) size += field_xy_.footprint();
    if ( field_partition_ ) size += field_partition_.footprint();
    if ( field_global_index_ ) size += field_global_index_.footprint();
    if ( field_remote_index_ ) size += field_remote_index_.footprint();
    if ( field_index_i_ ) size += field_index_i_.footprint();
    if ( field_index_j_ ) size += field_index_j_.footprint();
    return size;
}

}  // namespace detail

// ----------------------------------------------------------------------------

StructuredColumns::StructuredColumns() : FunctionSpace(), functionspace_( nullptr ) {}

StructuredColumns::StructuredColumns( const FunctionSpace& functionspace ) :
    FunctionSpace( functionspace ),
    functionspace_( dynamic_cast<const detail::StructuredColumns*>( get() ) ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const eckit::Configuration& config ) :
    FunctionSpace( new detail::StructuredColumns( grid, config ) ),
    functionspace_( dynamic_cast<const detail::StructuredColumns*>( get() ) ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Partitioner& partitioner,
                                      const eckit::Configuration& config ) :
    FunctionSpace( new detail::StructuredColumns( grid, partitioner, config ) ),
    functionspace_( dynamic_cast<const detail::StructuredColumns*>( get() ) ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const Vertical& vertical, const eckit::Configuration& config ) :
    FunctionSpace( new detail::StructuredColumns( grid, vertical, config ) ),
    functionspace_( dynamic_cast<const detail::StructuredColumns*>( get() ) ) {}

StructuredColumns::StructuredColumns( const Grid& grid, const Vertical& vertical, const grid::Partitioner& partitioner,
                                      const eckit::Configuration& config ) :
    FunctionSpace( new detail::StructuredColumns( grid, vertical, partitioner, config ) ),
    functionspace_( dynamic_cast<const detail::StructuredColumns*>( get() ) ) {}

void StructuredColumns::gather( const FieldSet& local, FieldSet& global ) const {
    functionspace_->gather( local, global );
}

void StructuredColumns::gather( const Field& local, Field& global ) const {
    functionspace_->gather( local, global );
}

void StructuredColumns::scatter( const FieldSet& global, FieldSet& local ) const {
    functionspace_->scatter( global, local );
}

void StructuredColumns::scatter( const Field& global, Field& local ) const {
    functionspace_->scatter( global, local );
}

std::string StructuredColumns::checksum( const FieldSet& fieldset ) const {
    return functionspace_->checksum( fieldset );
}

std::string StructuredColumns::checksum( const Field& field ) const {
    return functionspace_->checksum( field );
}

// ----------------------------------------------------------------------------

}  // namespace functionspace
}  // namespace atlas
