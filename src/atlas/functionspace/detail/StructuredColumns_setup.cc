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

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/parallel/mpi/mpi.h"
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
        idx_t r = static_cast<idx_t>( set.size() );
        set.emplace_back( i, j, r );
        return true;
    }

    idx_t size() const { return static_cast<idx_t>( set.size() ); }

    using const_iterator = decltype( set )::const_iterator;

    const_iterator begin() const { return set.begin(); }
    const_iterator end() const { return set.end(); }
};

}  // namespace


void StructuredColumns::setup( const grid::Distribution& distribution, const eckit::Configuration& config ) {
    ATLAS_TRACE( "Generating StructuredColumns..." );
    bool periodic_points = config.getInt( "periodic_points", false );
    if ( not( *grid_ ) ) {
        throw_Exception( "Grid is not a grid::Structured type", Here() );
    }
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
                    if ( periodic_points && i == grid_->nx( j ) - 1 ) {
                        ++i;
                    }

                    double x      = grid_->x( i, j );
                    double x_next = grid_->x( i + 1, j );
                    double x_prev = grid_->x( i - 1, j );
                    for ( idx_t jj = j - halo; jj <= j + halo; ++jj ) {
                        idx_t last = grid_->nx( compute_j( jj ) ) - 1;
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

            ATLAS_ASSERT( gridpoints.size() == owned );

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

            ATLAS_ASSERT( gridpoints.size() == owned + extra_halo );
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

        for ( const GridPoint& gp : gridpoints ) {
            xy( gp.r, XX ) = compute_x( gp.i, gp.j );
            if ( gp.j >= 0 && gp.j < grid_->ny() ) {
                xy( gp.r, YY ) = grid_->y( gp.j );
            }
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
            ghost( gp.r )   = 0;
        }

        for ( idx_t j = j_begin_halo_; j < j_begin_; ++j ) {
            for ( idx_t i = i_begin_halo_( j ); i < i_end_halo_( j ); ++i ) {
                ghost( index( i, j ) ) = 1;
            }
        }
        for ( idx_t j = j_begin_; j < j_end_; ++j ) {
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
