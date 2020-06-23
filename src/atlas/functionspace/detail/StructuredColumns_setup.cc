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
#include <numeric>
#include <sstream>
#include <string>

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace functionspace {
namespace detail {

namespace {

struct GridPoint {
public:
    idx_t i, j;
    idx_t r;

    GridPoint() = default;
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

    atlas::vector<GridPoint> gp_;

    void reserve( size_t size ) { gp_.reserve( size ); }
    void resize( size_t size ) { gp_.resize( size ); }
    void set( idx_t i, idx_t j, idx_t r ) {
        gp_[r].i = i;
        gp_[r].j = j;
        gp_[r].r = r;
    }

    idx_t size() const { return static_cast<idx_t>( gp_.size() ); }

    const GridPoint& operator[]( idx_t i ) const { return gp_[i]; }

    // unused:
    //using const_iterator = decltype( gp_ )::const_iterator;
    //const_iterator begin() const { return gp_.begin(); }
    //const_iterator end() const { return gp_.end(); }
};

}  // namespace


void StructuredColumns::setup( const grid::Distribution& distribution, const eckit::Configuration& config ) {
    config.get( "periodic_points", periodic_points_ );
    if ( not( *grid_ ) ) {
        throw_Exception( "Grid is not a grid::Structured type", Here() );
    }

    const double eps = 1.e-12;

    ny_                  = grid_->ny();
    north_pole_included_ = 90. - grid_->y( 0 ) == 0.;
    south_pole_included_ = 90. + grid_->y( ny_ - 1 ) == 0;

    distribution_ = distribution.type();

    nb_partitions_ = distribution.nb_partitions();

    int mpi_rank = int( mpi::rank() );
    int mpi_size = int( mpi::size() );

    j_begin_ = std::numeric_limits<idx_t>::max();
    j_end_   = std::numeric_limits<idx_t>::min();
    i_begin_.resize( grid_->ny(), std::numeric_limits<idx_t>::max() );
    i_end_.resize( grid_->ny(), std::numeric_limits<idx_t>::min() );
    idx_t owned( 0 );

    ATLAS_TRACE_SCOPE( "Compute bounds owned" ) {
        if ( mpi_size == 1 ) {
            j_begin_ = 0;
            j_end_   = grid_->ny();
            for ( idx_t j = j_begin_; j < j_end_; ++j ) {
                i_begin_[j] = 0;
                i_end_[j]   = grid_->nx( j );
            }
            owned = grid_->size();
        }
        else {
            size_t num_threads = atlas_omp_get_max_threads();
            if ( num_threads == 1 ) {
                idx_t c( 0 );
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
            }
            else {
                std::vector<idx_t> thread_reduce_j_begin( num_threads, std::numeric_limits<idx_t>::max() );
                std::vector<idx_t> thread_reduce_j_end( num_threads, std::numeric_limits<idx_t>::min() );
                std::vector<std::vector<idx_t>> thread_reduce_i_begin(
                    grid_->ny(), std::vector<idx_t>( num_threads, std::numeric_limits<idx_t>::max() ) );
                std::vector<std::vector<idx_t>> thread_reduce_i_end(
                    grid_->ny(), std::vector<idx_t>( num_threads, std::numeric_limits<idx_t>::min() ) );
                std::vector<idx_t> thread_reduce_owned( num_threads, 0 );
                atlas_omp_parallel {
                    const idx_t thread_num = atlas_omp_get_thread_num();
                    const idx_t begin =
                        static_cast<int>( size_t( thread_num ) * size_t( grid_->size() ) / size_t( num_threads ) );
                    const idx_t end =
                        static_cast<int>( size_t( thread_num + 1 ) * size_t( grid_->size() ) / size_t( num_threads ) );
                    idx_t thread_j_begin = 0;
                    std::vector<idx_t> thread_i_begin( grid_->ny() );
                    std::vector<idx_t> thread_i_end( grid_->ny() );
                    idx_t n = 0;
                    for ( idx_t j = 0; j < grid_->ny(); ++j ) {
                        if ( n + grid_->nx( j ) > begin ) {
                            thread_j_begin    = j;
                            thread_i_begin[j] = begin - n;
                            break;
                        }
                        n += grid_->nx( j );
                    }
                    idx_t thread_j_end{thread_j_begin};
                    for ( idx_t j = thread_j_begin; j < grid_->ny(); ++j ) {
                        idx_t i_end = end - n;
                        if ( j > thread_j_begin ) {
                            thread_i_begin[j] = 0;
                        }
                        if ( i_end > grid_->nx( j ) ) {
                            thread_i_end[j] = grid_->nx( j );
                            n += grid_->nx( j );
                        }
                        else {
                            thread_i_end[j] = i_end;
                            thread_j_end    = j + 1;
                            break;
                        }
                    }
                    auto& __j_begin = thread_reduce_j_begin[thread_num];
                    auto& __j_end   = thread_reduce_j_end[thread_num];
                    auto& __owned   = thread_reduce_owned[thread_num];
                    idx_t c         = begin;

                    for ( idx_t j = thread_j_begin; j < thread_j_end; ++j ) {
                        auto& __i_begin = thread_reduce_i_begin[j][thread_num];
                        auto& __i_end   = thread_reduce_i_end[j][thread_num];
                        bool j_in_partition{false};
                        for ( idx_t i = thread_i_begin[j]; i < thread_i_end[j]; ++i, ++c ) {
                            if ( distribution.partition( c ) == mpi_rank ) {
                                j_in_partition = true;
                                __i_begin      = std::min<idx_t>( __i_begin, i );
                                __i_end        = std::max<idx_t>( __i_end, i + 1 );
                                ++__owned;
                            }
                        }
                        if ( j_in_partition ) {
                            __j_begin = std::min<idx_t>( __j_begin, j );
                            __j_end   = std::max<idx_t>( __j_end, j + 1 );
                        }
                    }
                    ATLAS_ASSERT( c == end );
                }
                owned    = std::accumulate( thread_reduce_owned.begin(), thread_reduce_owned.end(), 0 );
                j_begin_ = *std::min_element( thread_reduce_j_begin.begin(), thread_reduce_j_begin.end() );
                j_end_   = *std::max_element( thread_reduce_j_end.begin(), thread_reduce_j_end.end() );
                atlas_omp_parallel_for( idx_t j = j_begin_; j < j_end_; ++j ) {
                    i_begin_[j] = *std::min_element( thread_reduce_i_begin[j].begin(), thread_reduce_i_begin[j].end() );
                    i_end_[j]   = *std::max_element( thread_reduce_i_end[j].begin(), thread_reduce_i_end[j].end() );
                }
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

    auto compute_i_fast = []( idx_t i, const idx_t nx ) -> idx_t {
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
        if ( j < 0 ) {
            j = ( grid_->y( 0 ) == 90. ) ? -j : -j - 1;
        }
        else if ( j >= grid_->ny() ) {
            idx_t jlast = grid_->ny() - 1;
            j           = ( grid_->y( jlast ) == -90. ) ? jlast - 1 - ( j - grid_->ny() ) : jlast - ( j - grid_->ny() );
        }
        if ( j < 0 or j >= grid_->ny() ) {
            j = compute_j( j );
        }
        return j;
    };

    auto compute_x = [this, &compute_i, &compute_j]( idx_t i, idx_t j ) -> double {
        const idx_t jj = compute_j( j );
        const idx_t ii = compute_i( i, jj );  // guaranteed between 0 and nx(jj)
        const idx_t nx = grid_->nx( jj );
        const double a = ( ii - i ) / nx;
        const double x = grid_->x( ii, jj ) - a * grid_->x( nx, jj );
        return x;
    };

    auto compute_x_fast = [this, &compute_i_fast]( idx_t i, idx_t jj, idx_t nx ) -> double {
        const idx_t ii = compute_i_fast( i, nx );  // guaranteed between 0 and nx(jj)
        const double a = ( ii - i ) / nx;
        const double x = grid_->x( ii, jj ) - a * ( grid_->x( nx, jj ) - grid_->x( 0, jj ) );
        return x;
    };

    auto compute_i_less_equal_x = [this, &eps]( const double& x, idx_t j ) -> idx_t {
        const double dx = grid_->dx( j );
        idx_t i         = idx_t( std::floor( ( x + eps - grid_->xmin( j ) ) / dx ) );
        return i;
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
            ATLAS_ASSERT( grid_->nx( jj ) % 2 == 0 );  // assert even number of points
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
            ATLAS_ASSERT( grid_->nx( jj ) % 2 == 0 );  // assert even number of points
            ii = ( ii < grid_->nx( jj ) / 2 ) ? ii + grid_->nx( jj ) / 2
                                              : ( ii >= grid_->nx( jj ) / 2 ) ? ii - grid_->nx( jj ) / 2 : ii;
        }
        p = distribution.partition( global_offsets[jj] + ii );
        return p;
    };

    GridPointSet gridpoints;

    ATLAS_TRACE_SCOPE( "Compute mapping" ) {
        idx_t imin = std::numeric_limits<idx_t>::max();
        idx_t imax = -std::numeric_limits<idx_t>::max();
        idx_t jmin = std::numeric_limits<idx_t>::max();
        idx_t jmax = -std::numeric_limits<idx_t>::max();

        ATLAS_TRACE_SCOPE( "Compute bounds halo" ) {
            for ( idx_t j = j_begin_halo_; j < j_end_halo_; ++j ) {
                i_begin_halo_( j ) = imin;
                i_end_halo_( j )   = imax;
            }

            // Following cannot be multithreaded in current form due to race-conditions related to index jj
            for ( idx_t j = j_begin_; j < j_end_; ++j ) {
                for ( idx_t i : {i_begin_[j], i_end_[j] - 1} ) {
                    // Following line only, increases periodic halo on the east side by 1
                    if ( periodic_points_ && i == grid_->nx( j ) - 1 ) {
                        ++i;
                    }

                    double x = grid_->x( i, j );

                    double x_next = grid_->x( i + 1, j );
                    double x_prev = grid_->x( i - 1, j );
                    for ( idx_t jj = j - halo; jj <= j + halo; ++jj ) {
                        idx_t jjj    = compute_j( jj );
                        idx_t nx_jjj = grid_->nx( jjj );
                        idx_t last   = grid_->nx( jjj ) - 1;
                        if ( i == grid_->nx( j ) ) {
                            ++last;
                        }

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

                        // idx_t ii = -halo;
                        // while ( compute_x_fast( ii, jjj, nx_jjj ) < x - eps ) {
                        //     ii++;
                        // }
                        // Question: is following implementation reproducible with above original while loop?
                        idx_t ii = compute_i_less_equal_x( x, jjj );

                        // ATLAS-186 workaround
                        // This while should not have to be there, but is here because of
                        // the MatchingMeshDomainDecomposition algorithm. that may lead to points
                        // left of the point ii.
                        while ( compute_x_fast( ii - 1, jjj, nx_jjj ) > x_prev + eps ) {
                            --ii;
                        }

                        idx_t i_minus_halo = ii - halo;

                        // Compute iii as index less-equal of x_next
                        //
                        //               x(i,j) x_next(i,j)
                        //   |-----|-----|-+---|-+---|-----|
                        // ii-halo            iii       iii+halo
                        //
                        idx_t iii = ii;
                        while ( compute_x_fast( iii + 1, jjj, nx_jjj ) < x_next - eps ) {
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


        ATLAS_TRACE_SCOPE( "Assemble gridpoints" ) {
            gridpoints.reserve( owned + extra_halo );
            gridpoints.resize( owned );

            if ( atlas_omp_get_max_threads() == 1 ) {
                idx_t r = 0;
                for ( idx_t j = j_begin_; j < j_end_; ++j ) {
                    for ( idx_t i = i_begin_[j]; i < i_end_[j]; ++i, ++r ) {
                        gridpoints.set( i, j, r );
                    }
                }
                ATLAS_ASSERT( r == owned );
            }
            else {
                atlas_omp_parallel {
                    const size_t num_threads = atlas_omp_get_num_threads();
                    const size_t thread_num  = atlas_omp_get_thread_num();

                    auto chunksize = size_t( owned ) / num_threads;
                    size_t begin   = chunksize * thread_num;
                    size_t end     = ( thread_num == num_threads - 1 ) ? owned : begin + chunksize;

                    idx_t thread_j_begin = 0;
                    std::vector<idx_t> thread_i_begin( grid_->ny() );
                    std::vector<idx_t> thread_i_end( grid_->ny() );
                    size_t n = 0;
                    for ( idx_t j = j_begin_; j < j_end_; ++j ) {
                        size_t j_size = static_cast<size_t>( i_end_[j] - i_begin_[j] );
                        if ( n + j_size > begin ) {
                            thread_j_begin    = j;
                            thread_i_begin[j] = i_begin_[j] + static_cast<idx_t>( begin - n );
                            break;
                        }
                        n += j_size;
                    }
                    idx_t thread_j_end{thread_j_begin};
                    for ( idx_t j = thread_j_begin; j < j_end_; ++j ) {
                        if ( j > thread_j_begin ) {
                            thread_i_begin[j] = i_begin_[j];
                        }
                        idx_t j_size = ( i_end_[j] - i_begin_[j] );

                        idx_t remaining = static_cast<idx_t>( end - n );
                        if ( j_size <= remaining ) {
                            thread_i_end[j] = i_end_[j];
                            n += j_size;
                            if ( n == end ) {
                                goto stop;
                            }
                        }
                        else {
                            thread_i_end[j] = i_begin_[j] + remaining;
                            goto stop;
                        }

                        continue;

                    stop:
                        thread_j_end = j + 1;
                        break;
                    }

                    idx_t r = idx_t( begin );
                    for ( idx_t j = thread_j_begin; j < thread_j_end; ++j ) {
                        for ( idx_t i = thread_i_begin[j]; i < thread_i_end[j]; ++i, ++r ) {
                            gridpoints.set( i, j, r );
                        }
                    }
                    if ( r != idx_t( end ) ) {
                        ATLAS_DEBUG_VAR( thread_num );
                        ATLAS_DEBUG_VAR( begin );
                        ATLAS_DEBUG_VAR( end );
                        ATLAS_DEBUG_VAR( r );
                    }
                    ATLAS_ASSERT( r == idx_t( end ) );
                }
            }

            ATLAS_ASSERT( gridpoints.size() == owned );

            gridpoints.resize( owned + extra_halo );
            idx_t r = owned;
            for ( idx_t j = j_begin_halo_; j < j_begin_; ++j ) {
                for ( idx_t i = i_begin_halo_( j ); i < i_end_halo_( j ); ++i, ++r ) {
                    gridpoints.set( i, j, r );
                }
            }
            for ( idx_t j = j_begin_; j < j_end_; ++j ) {
                for ( idx_t i = i_begin_halo_( j ); i < i_begin_[j]; ++i, ++r ) {
                    gridpoints.set( i, j, r );
                }
                for ( idx_t i = i_end_[j]; i < i_end_halo_( j ); ++i, ++r ) {
                    gridpoints.set( i, j, r );
                }
            }
            for ( idx_t j = j_end_; j < j_end_halo_; ++j ) {
                for ( idx_t i = i_begin_halo_( j ); i < i_end_halo_( j ); ++i, ++r ) {
                    gridpoints.set( i, j, r );
                }
            }

            ATLAS_ASSERT( gridpoints.size() == owned + extra_halo );
        }

        ATLAS_TRACE_SCOPE( "Fill in ij2gp " ) {
            ij2gp_.resize( {imin, imax}, {jmin, jmax} );

            atlas_omp_parallel_for( idx_t n = 0; n < gridpoints.size(); ++n ) {
                const GridPoint& gp = gridpoints[n];
                ij2gp_.set( gp.i, gp.j, gp.r );
            }
        }
    }

    ATLAS_TRACE_SCOPE( "Create fields" ) {
        size_halo_          = gridpoints.size();
        field_partition_    = Field( "partition", array::make_datatype<int>(), array::make_shape( size_halo_ ) );
        field_ghost_        = Field( "ghost", array::make_datatype<int>(), array::make_shape( size_halo_ ) );
        field_global_index_ = Field( "glb_idx", array::make_datatype<gidx_t>(), array::make_shape( size_halo_ ) );
        field_index_i_      = Field( "index_i", array::make_datatype<idx_t>(), array::make_shape( size_halo_ ) );
        field_index_j_      = Field( "index_j", array::make_datatype<idx_t>(), array::make_shape( size_halo_ ) );
        field_xy_           = Field( "xy", array::make_datatype<double>(), array::make_shape( size_halo_, 2 ) );

        auto xy         = array::make_view<double, 2>( field_xy_ );
        auto part       = array::make_view<int, 1>( field_partition_ );
        auto ghost      = array::make_view<int, 1>( field_ghost_ );
        auto global_idx = array::make_view<gidx_t, 1>( field_global_index_ );
        auto index_i    = array::make_indexview<idx_t, 1>( field_index_i_ );
        auto index_j    = array::make_indexview<idx_t, 1>( field_index_j_ );

        atlas_omp_parallel_for( idx_t n = 0; n < gridpoints.size(); ++n ) {
            const GridPoint& gp = gridpoints[n];
            if ( gp.j >= 0 && gp.j < grid_->ny() ) {
                xy( gp.r, XX ) = grid_->x( gp.i, gp.j );
                xy( gp.r, YY ) = grid_->y( gp.j );
            }
            else {
                xy( gp.r, XX ) = compute_x( gp.i, gp.j );
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
            ghost( gp.r )   = 0;
        }

        // Following short loops are not parallelized with

        for ( idx_t j = j_begin_halo_; j < j_begin_; ++j ) {
            for ( idx_t i = i_begin_halo_( j ); i < i_end_halo_( j ); ++i ) {
                ghost( index( i, j ) ) = 1;
            }
        }
        atlas_omp_parallel_for( idx_t j = j_begin_; j < j_end_; ++j ) {
            for ( idx_t i = i_begin_halo_( j ); i < i_begin_[j]; ++i ) {
                ghost( index( i, j ) ) = 1;
            }
            for ( idx_t i = i_end_[j]; i < i_end_halo_( j ); ++i ) {
                ghost( index( i, j ) ) = 1;
            }
        }
        for ( idx_t j = j_end_; j < j_end_halo_; ++j ) {
            for ( idx_t i = i_begin_halo_( j ); i < i_end_halo_( j ); ++i ) {
                ghost( index( i, j ) ) = 1;
            }
        }
    }
}

}  // namespace detail

// ----------------------------------------------------------------------------

}  // namespace functionspace
}  // namespace atlas
