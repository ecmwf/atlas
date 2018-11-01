/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include <cmath>
#include <limits>

#include "atlas/interpolation/method/structured/Bicubic.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Point3.h"
#include "eckit/log/Plural.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/log/Seconds.h"
#include "eckit/mpi/Comm.h"

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/Buffer.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"


#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace {

MethodBuilder<Bicubic> __builder1( "structured-bicubic" );
MethodBuilder<Bicubic> __builder2( "bicubic" );

}  // namespace

void Bicubic::setup( const FunctionSpace& source, const FunctionSpace& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::CubicStructured2D::setup()" );

    source_ = source;
    target_ = target;

    if ( functionspace::NodeColumns tgt = target ) {
        target_xy_    = tgt.mesh().nodes().lonlat();
        target_ghost_ = tgt.mesh().nodes().ghost();
    }
    else if ( functionspace::PointCloud tgt = target ) {
        target_xy_    = tgt.lonlat();
        target_ghost_ = tgt.ghost();
    }
    else {
        NOTIMP;
    }

    setup( source );
}

void Bicubic::print( std::ostream& out ) const {
    ASSERT( not matrix_.empty() );

    functionspace::NodeColumns src( source_ );
    functionspace::NodeColumns tgt( target_ );
    if ( not tgt ) NOTIMP;
    auto gidx_src = array::make_view<gidx_t, 1>( src.nodes().global_index() );

    ASSERT( tgt.nodes().size() == idx_t( matrix_.rows() ) );


    auto field_stencil_points_loc  = tgt.createField<gidx_t>( option::variables( 16 ) );
    auto field_stencil_weigths_loc = tgt.createField<double>( option::variables( 16 ) );
    auto field_stencil_size_loc    = tgt.createField<int>();

    auto stencil_points_loc  = array::make_view<gidx_t, 2>( field_stencil_points_loc );
    auto stencil_weights_loc = array::make_view<double, 2>( field_stencil_weigths_loc );
    auto stencil_size_loc    = array::make_view<int, 1>( field_stencil_size_loc );
    stencil_size_loc.assign( 0 );

    for ( Matrix::const_iterator it = matrix_.begin(); it != matrix_.end(); ++it ) {
        idx_t p                     = idx_t( it.row() );
        idx_t& i                    = stencil_size_loc( p );
        stencil_points_loc( p, i )  = gidx_src( it.col() );
        stencil_weights_loc( p, i ) = *it;
        ++i;
    }


    gidx_t global_size = tgt.gather().glb_dof();

    auto field_stencil_points_glb  = tgt.createField<gidx_t>( option::variables( 16 ) | option::global( 0 ) );
    auto field_stencil_weights_glb = tgt.createField<double>( option::variables( 16 ) | option::global( 0 ) );
    auto field_stencil_size_glb    = tgt.createField<int>( option::global( 0 ) );


    auto stencil_points_glb  = array::make_view<gidx_t, 2>( field_stencil_points_glb );
    auto stencil_weights_glb = array::make_view<double, 2>( field_stencil_weights_glb );
    auto stencil_size_glb    = array::make_view<int, 1>( field_stencil_size_glb );

    tgt.gather().gather( stencil_size_loc, stencil_size_glb );
    tgt.gather().gather( stencil_points_loc, stencil_points_glb );
    tgt.gather().gather( stencil_weights_loc, stencil_weights_glb );

    if ( mpi::comm().rank() == 0 ) {
        int precision = std::numeric_limits<double>::max_digits10;
        for ( idx_t i = 0; i < global_size; ++i ) {
            out << std::setw( 10 ) << i + 1 << " : ";
            for ( idx_t j = 0; j < stencil_size_glb( i ); ++j ) {
                out << std::setw( 10 ) << stencil_points_glb( i, j );
            }
            for ( idx_t j = stencil_size_glb( i ); j < 16; ++j ) {
                out << "          ";
            }
            for ( idx_t j = 0; j < stencil_size_glb( i ); ++j ) {
                out << std::setw( precision + 5 ) << std::left << std::setprecision( precision )
                    << stencil_weights_glb( i, j );
            }
            out << std::endl;
        }
    }
}

void Bicubic::setup( const FunctionSpace& source ) {
    src_ = source;
    ASSERT( src_ );
    ASSERT( src_.halo() >= 2 );

    compute_horizontal_stencil_ = ComputeHorizontalStencil( src_.grid(), stencil_width() );

    if ( not matrix_free_ ) {
        idx_t inp_npts = src_.size();
        idx_t out_npts = target_xy_.shape( 0 );

        auto ghost = array::make_view<int, 1>( target_ghost_ );
        auto xy    = array::make_view<double, 2>( target_xy_ );

        auto triplets = reserve_triplets( out_npts );

        WorkSpace workspace;
        ATLAS_TRACE_SCOPE( "Computing interpolation matrix" ) {
            eckit::ProgressTimer progress( "Computing interpolation weights", out_npts, "point", double( 5 ),
                                           Log::debug() );
            for ( idx_t n = 0; n < out_npts; ++n, ++progress ) {
                PointXY p{xy( n, XX ), xy( n, YY )};
                while ( p.x() < 0. ) {
                    p.x() += 360.;
                }
                while ( p.x() >= 360. ) {
                    p.x() -= 360.;
                }
                if ( not ghost( n ) ) { insert_triplets( n, p, triplets, workspace ); }
            }
        }

        // fill sparse matrix and return
        Matrix A( out_npts, inp_npts, triplets );
        matrix_.swap( A );
    }
    else {
        NOTIMP;
    }
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
