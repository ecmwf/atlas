/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"

namespace atlas {
namespace meshgenerator {
namespace detail {
namespace cubedsphere {


// -----------------------------------------------------------------------------
// Projection cast
// -----------------------------------------------------------------------------

const CubedSphereProjectionBase* castProjection( const ProjectionImpl* projectionPtr ) {
    const auto* const csProjectionPtr = dynamic_cast<const CubedSphereProjectionBase*>( projectionPtr );

    if ( !csProjectionPtr ) {
        throw_Exception( "Cannot cast ProjectionImpl* to CubedSphereProjectionBase*.", Here() );
    }
    return csProjectionPtr;
}

// -----------------------------------------------------------------------------
// Jacobian2 class
// -----------------------------------------------------------------------------

Jacobian2::Jacobian2( double df0_by_dx0, double df0_by_dx1, double df1_by_dx0, double df1_by_dx1 ) :
    df0_by_dx0_( df0_by_dx0 ), df0_by_dx1_( df0_by_dx1 ), df1_by_dx0_( df1_by_dx0 ), df1_by_dx1_( df1_by_dx1 ) {}

Jacobian2::Jacobian2( const Point2& f00, const Point2& f10, const Point2& f01, double dx0, double dx1 ) :
    Jacobian2( ( f10[0] - f00[0] ) / dx0, ( f01[0] - f00[0] ) / dx1, ( f10[1] - f00[1] ) / dx0,
               ( f01[1] - f00[1] ) / dx1 ) {}

Jacobian2::Jacobian2( const Point2& f00, const Point2& f10, const Point2& f01 ) : Jacobian2( f00, f10, f01, 1., 1. ) {}

double Jacobian2::det() const {
    return df0_by_dx0_ * df1_by_dx1_ - df0_by_dx1_ * df1_by_dx0_;
}

Jacobian2 Jacobian2::operator*( double a ) const {
    return Jacobian2( df0_by_dx0_ * a, df0_by_dx1_ * a, df1_by_dx0_ * a, df1_by_dx1_ * a );
}

Point2 Jacobian2::operator*( const Point2& dx ) const {
    return Point2( dx[0] * df0_by_dx0_ + df0_by_dx1_ * dx[1], dx[0] * df1_by_dx0_ + df1_by_dx1_ * dx[1] );
}

Jacobian2 Jacobian2::operator*( const Jacobian2& Jb ) const {
    return Jacobian2( df0_by_dx0_ * Jb.df0_by_dx0_ + df0_by_dx1_ * Jb.df1_by_dx0_,
                      df0_by_dx0_ * Jb.df0_by_dx1_ + df0_by_dx1_ * Jb.df1_by_dx1_,
                      df1_by_dx0_ * Jb.df0_by_dx0_ + df1_by_dx1_ * Jb.df1_by_dx0_,
                      df1_by_dx0_ * Jb.df0_by_dx1_ + df1_by_dx1_ * Jb.df1_by_dx1_ );
}

Jacobian2 Jacobian2::inverse() const {
    return Jacobian2( df1_by_dx1_, -df0_by_dx1_, -df1_by_dx0_, df0_by_dx0_ ) * ( 1. / det() );
}

Jacobian2 Jacobian2::sign() const {
    const double smallNumber = det() * std::numeric_limits<double>::epsilon();
    const auto signValue     = [&]( double number ) -> double {
        return std::abs( number ) < smallNumber ? 0. : number < 0. ? -1. : 1.;
    };
    return Jacobian2( signValue( df0_by_dx0_ ), signValue( df0_by_dx1_ ), signValue( df1_by_dx0_ ),
                      signValue( df1_by_dx1_ ) );
}


// -----------------------------------------------------------------------------
// NeighbourJacobian class
// -----------------------------------------------------------------------------

NeighbourJacobian::NeighbourJacobian( const CubedSphereGrid& csGrid ) {
    // N must be greater than 2.
    if ( csGrid.N() < 2 ) {
        throw_Exception( "Jacobians can only be calculated for N > 1 .", Here() );
    }

    // Assumes cell centre staggering.
    if ( csGrid.stagger() != "C" ) {
        throw_Exception( "NeighbourJacobian class only works for cell-centre grid", Here() );
    }

    // Get projection.
    csProjection_ = castProjection( csGrid.projection().get() );

    // Get tiles.
    const auto& csTiles = csProjection_->getCubedSphereTiles();

    // Get grid size.
    N_ = csGrid.N();


    // Get cell width.
    const double cellWidth = 90. / N_;

    // Get xy of points (i = 0, j = 0), (i = 1, j = 0) and (i = 0, j = 1) on tiles.
    std::array<PointXY, 6> xy00;
    std::array<PointXY, 6> xy10;
    std::array<PointXY, 6> xy01;

    // Loop over grid points.
    auto tijIt = csGrid.tij().begin();
    for ( const PointXY& xy : csGrid.xy() ) {
        const auto t  = static_cast<size_t>( ( *tijIt ).t() );
        const idx_t i = ( *tijIt ).i();
        const idx_t j = ( *tijIt ).j();


        if ( i == 0 && j == 0 ) {
            xy00[t] = xy;
        }
        else if ( i == 1 && j == 0 ) {
            xy10[t] = xy;
        }
        else if ( i == 0 && j == 1 ) {
            xy01[t] = xy;
        }
        ++tijIt;
    }

    for ( size_t t = 0; t < 6; ++t ) {
        // Calculate tile Jacobians.
        dxy_by_dij_[t] = Jacobian2( xy00[t], xy10[t], xy01[t] );

        // Rescale by cell width (gains an extra couple of decimal places of precision).
        dxy_by_dij_[t] = dxy_by_dij_[t].sign() * cellWidth;

        // Get inverse.
        dij_by_dxy_[t] = dxy_by_dij_[t].inverse();

        // Set xy00. Grid point needs moving to (i = 0, j = 0).
        xy00_[t] = xy00[t] + dxy_by_dij_[t] * PointIJ( -0.5, -0.5 );

        // Get other three corners so we can work out xy min/max.
        const PointXY xyN0 = xy00_[t] + dxy_by_dij_[t] * PointIJ( N_, 0 );
        const PointXY xyNN = xy00_[t] + dxy_by_dij_[t] * PointIJ( N_, N_ );
        const PointXY xy0N = xy00_[t] + dxy_by_dij_[t] * PointIJ( 0, N_ );

        // Get xy min/max.
        std::tie( xyMin_[t].x(), xyMax_[t].x() ) = std::minmax( {xy00_[t].x(), xyN0.x(), xyNN.x(), xy0N.x()} );
        std::tie( xyMin_[t].y(), xyMax_[t].y() ) = std::minmax( {xy00_[t].y(), xyN0.y(), xyNN.y(), xy0N.y()} );

        // Round to nearest degree.
        xyMin_[t].x() = std::round( xyMin_[t].x() );
        xyMax_[t].x() = std::round( xyMax_[t].x() );
        xyMin_[t].y() = std::round( xyMin_[t].y() );
        xyMax_[t].y() = std::round( xyMax_[t].y() );

        // Neighbour assignment lambda.
        const auto neighbourAssignment = [&]( TileEdge::k k ) -> void {
            // Shift points in to neighbouring tiles.
            PointIJ ijDisplacement;
            switch ( k ) {
                case TileEdge::LEFT: {
                    ijDisplacement = PointIJ( -2, 0 );
                    break;
                }
                case TileEdge::BOTTOM: {
                    ijDisplacement = PointIJ( 0, -2 );
                    break;
                }
                case TileEdge::RIGHT: {
                    ijDisplacement = PointIJ( N_, 0 );
                    break;
                }
                case TileEdge::TOP: {
                    ijDisplacement = PointIJ( 0, N_ );
                    break;
                }
                case TileEdge::UNDEFINED: {
                    throw_Exception( "Undefined tile edge.", Here() );
                }
            }

            // Convert displacement from ij to xy.
            const PointXY xyDisplacement = dxy_by_dij_[t] * ijDisplacement;

            // Get neighbour xy points in xy space local to tile.
            const PointXY xy00Local = xy00[t] + xyDisplacement;
            const PointXY xy10Local = xy10[t] + xyDisplacement;
            const PointXY xy01Local = xy01[t] + xyDisplacement;

            // Convert from local xy to global xy.
            const PointXY xy00Global = csTiles.tileCubePeriodicity( xy00Local, static_cast<idx_t>( t ) );
            const PointXY xy10Global = csTiles.tileCubePeriodicity( xy10Local, static_cast<idx_t>( t ) );
            const PointXY xy01Global = csTiles.tileCubePeriodicity( xy01Local, static_cast<idx_t>( t ) );

            // Get neighbour tile ID.
            neighbours_[t].t_[k] = csTiles.indexFromXY( xy00Global.data() );

            // Set Jacobian of global xy with respect to local ij.
            auto dxyGlobal_by_dij = Jacobian2( xy00Global, xy10Global, xy01Global );

            // Rescale by cell width (gains an extra couple of decimal places of precision).
            dxyGlobal_by_dij = dxyGlobal_by_dij.sign() * cellWidth;

            // Chain rule to get Jacobian with respect to local xy.
            neighbours_[t].dxyGlobal_by_dxyLocal_[k] = dxyGlobal_by_dij * dij_by_dxy_[t];

            // Set local xy00
            neighbours_[t].xy00Local_[k] = xy00Local;

            // Set global xy00
            neighbours_[t].xy00Global_[k] = xy00Global;
        };

        // Assign neighbours (good job we put it all in a lambda!).
        neighbourAssignment( TileEdge::LEFT );
        neighbourAssignment( TileEdge::BOTTOM );
        neighbourAssignment( TileEdge::RIGHT );
        neighbourAssignment( TileEdge::TOP );
    }
}

PointXY NeighbourJacobian::xy( const PointIJ& ij, idx_t t ) const {
    // Get jacobian.
    const Jacobian2& jac = dxy_by_dij_[static_cast<size_t>( t )];
    const PointXY& xy00  = xy00_[static_cast<size_t>( t )];

    // Return ij
    return xy00 + jac * ij;
}

PointIJ NeighbourJacobian::ij( const PointXY& xy, idx_t t ) const {
    // Get jacobian.
    const Jacobian2& jac = dij_by_dxy_[static_cast<size_t>( t )];
    const PointXY& xy00  = xy00_[static_cast<size_t>( t )];

    // Return ij
    return jac * ( xy - xy00 );
}

PointTXY NeighbourJacobian::xyLocalToGlobal( const PointXY& xyLocal, idx_t tLocal ) const {
    // The tileCubePeriodicity method fails when extrapolating along an unowned
    // tile edge. This method explicitly places an xy point on to a neighbouring
    // tile to avoid this. Once the correct xy position has been found,
    // tileCubePeriodicity will correcty find the "owned" xy position of a point
    // on an unowned tile edge.

    // Declare result.
    PointXY xyGlobal;
    idx_t tGlobal;

    // Get ij.
    const PointIJ ijLocal = ij( xyLocal, tLocal );

    // Get tiles.
    const auto& csTiles = csProjection_->getCubedSphereTiles();

    if ( ijInterior( ijLocal ) ) {
        // We're within the tile boundary (possibly on an edge).

        // Return local values if not on edge.
        if ( !ijEdge( ijLocal ) ) {
            return PointTXY( tLocal, xyLocal );
        }

        // We're on an edge. Will need to check with Tiles class.
        xyGlobal = xyLocal;
        tGlobal  = tLocal;
    }
    else {
        // We're outside the tile boundary.
        // Figure out which tile xy is on.
        TileEdge::k k;
        if ( ijLocal.iNode() < 0 ) {
            k = TileEdge::LEFT;
        }
        else if ( ijLocal.jNode() < 0 ) {
            k = TileEdge::BOTTOM;
        }
        else if ( ijLocal.iNode() > N_ ) {
            k = TileEdge::RIGHT;
        }
        else if ( ijLocal.jNode() > N_ ) {
            k = TileEdge::TOP;
        }
        else {
            throw_Exception( "Cannot determine neighbour tile.", Here() );
        }

        // Get reference points and jacobian.
        const PointXY& xy00Local_  = neighbours_[static_cast<size_t>( tLocal )].xy00Local_[k];
        const PointXY& xy00Global_ = neighbours_[static_cast<size_t>( tLocal )].xy00Global_[k];
        const Jacobian2& jac       = neighbours_[static_cast<size_t>( tLocal )].dxyGlobal_by_dxyLocal_[k];

        // Get t.
        tGlobal = neighbours_[static_cast<size_t>( tLocal )].t_[k];

        // Calculate global xy.
        xyGlobal = xy00Global_ + jac * ( xyLocal - xy00Local_ );
    }

    // Need to be very careful with floating point comparisons used in projection
    // class. Move point on to edge if it is very close.
    xyGlobal = snapToEdge( xyGlobal, tGlobal );

    // Correct for edge-ownership rules.
    xyGlobal = csTiles.tileCubePeriodicity( xyGlobal, tGlobal );
    tGlobal  = csTiles.indexFromXY( xyGlobal.data() );

    return PointTXY( tGlobal, xyGlobal );
}

PointTIJ NeighbourJacobian::ijLocalToGlobal( const PointIJ& ijLocal, idx_t tLocal ) const {
    // Use xyLocalToGlobal method to take care of this.

    // Get global xyt.
    PointTXY txyGlobal = xyLocalToGlobal( xy( ijLocal, tLocal ), tLocal );

    // convert to ijt
    return PointTIJ( txyGlobal.t(), ij( txyGlobal.xy(), txyGlobal.t() ) );
}

bool NeighbourJacobian::ijInterior( const PointIJ& ij ) const {
    return ij.iNode() >= 0 && ij.iNode() <= N_ && ij.jNode() >= 0 && ij.jNode() <= N_;
}

bool NeighbourJacobian::ijEdge( const PointIJ& ij ) const {
    return ijInterior( ij ) && ( ij.iNode() == 0 || ij.iNode() == N_ || ij.jNode() == 0 || ij.jNode() == N_ );
}

bool NeighbourJacobian::ijCross( const PointIJ& ij ) const {
    const bool inCorner = ( ij.iNode() < 0 && ij.jNode() < 0 ) ||    // bottom-left corner.
                          ( ij.iNode() > N_ && ij.jNode() < 0 ) ||   // bottom-right corner.
                          ( ij.iNode() > N_ && ij.jNode() > N_ ) ||  // top-right corner.
                          ( ij.iNode() < 0 && ij.jNode() > N_ );     // top-left corner.
    return !inCorner;
}

PointXY NeighbourJacobian::snapToEdge( const PointXY& xy, idx_t t ) const {
    const auto nudgeValue = []( double a, double b ) -> double {
        // Set tolerance to machine epsilon * 360 degrees.
        constexpr double tol = 360. * std::numeric_limits<double>::epsilon();

        // If a is nearly equal to b, return b. Otherwise return a.
        return std::abs( a - b ) <= tol ? b : a;
    };

    // If point is near edge, place it exactly on edge.
    const PointXY& xyMin = xyMin_[static_cast<size_t>( t )];
    const PointXY& xyMax = xyMax_[static_cast<size_t>( t )];
    return PointXY( nudgeValue( nudgeValue( xy.x(), xyMin.x() ), xyMax.x() ),
                    nudgeValue( nudgeValue( xy.y(), xyMin.y() ), xyMax.y() ) );
}

}  // namespace cubedsphere
}  // namespace detail
}  // namespace meshgenerator
}  // namespace atlas
