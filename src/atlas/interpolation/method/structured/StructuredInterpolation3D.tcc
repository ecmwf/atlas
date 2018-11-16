/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#pragma once

#include "StructuredInterpolation3D.h"


#include "eckit/exception/Exceptions.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {


template <typename Kernel>
StructuredInterpolation3D<Kernel>::StructuredInterpolation3D( const Method::Config& config ) :
    Method( config ),
    matrix_free_{false} {
    config.get( "matrix_free", matrix_free_ );

    if ( not matrix_free_ ) {
        throw eckit::NotImplemented( "Matrix-free StructuredInterpolation3D not implemented", Here() );
    }
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::setup( const Grid& source, const Grid& target ) {
    if ( mpi::comm().size() > 1 ) { NOTIMP; }


    ASSERT( grid::StructuredGrid( source ) );
    FunctionSpace source_fs = functionspace::StructuredColumns( source, option::halo( kernel_->stencil_halo() ) );
    FunctionSpace target_fs = functionspace::PointCloud( target );

    setup( source_fs, target_fs );
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::setup( const FunctionSpace& source, const FunctionSpace& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::StructuredInterpolation::setup()" );

    source_ = source;
    target_ = target;

    if ( functionspace::PointCloud tgt = target ) {
        target_lonlat_   = tgt.lonlat();
        target_vertical_ = tgt.vertical();
        target_ghost_    = tgt.ghost();
    }
    else {
        NOTIMP;
    }

    setup( source );
}

template <typename Kernel>
void StructuredInterpolation3D<Kernel>::setup( const FunctionSpace& source, const Field& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::StructuredInterpolation::setup()" );

    source_ = source;

    if ( target.functionspace() ) { target_ = target.functionspace(); }
    ASSERT( target.levels() );

    target_3d_ = target;

    setup( source );
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::print( std::ostream& ) const {
    NOTIMP;
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::setup( const FunctionSpace& source ) {
    kernel_.reset( new Kernel( source ) );

    if ( not matrix_free_ ) {
        idx_t inp_npts = source.size();
        idx_t out_npts = target_lonlat_.shape( 0 );

        auto ghost    = array::make_view<int, 1>( target_ghost_ );
        auto lonlat   = array::make_view<double, 2>( target_lonlat_ );
        auto vertical = array::make_view<double, 1>( target_vertical_ );

        auto triplets = kernel_->allocate_triplets( out_npts );

        constexpr NormaliseLongitude normalise;
        ATLAS_TRACE_SCOPE( "Precomputing interpolation matrix" ) {
            atlas_omp_parallel {
                typename Kernel::WorkSpace workspace;
                atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                    if ( not ghost( n ) ) {
                        PointLonLat p{lonlat( n, LON ), lonlat( n, LAT )};
                        normalise( p );
                        kernel_->insert_triplets( n, p, vertical( n ), triplets, workspace );
                    }
                }
            }
            // fill sparse matrix and return
            Matrix A( out_npts, inp_npts, triplets );
            matrix_.swap( A );
        }
    }
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::execute( const Field& src_field, Field& tgt_field ) const {
    FieldSet tgt( tgt_field );
    execute( FieldSet( src_field ), tgt );
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::execute( const FieldSet& src_fields, FieldSet& tgt_fields ) const {
    if ( not matrix_free_ ) {
        Method::execute( src_fields, tgt_fields );
        return;
    }

    ATLAS_TRACE( "atlas::interpolation::method::StructuredInterpolation::execute()" );

    const idx_t N = src_fields.size();
    ASSERT( N == tgt_fields.size() );

    if ( N == 0 ) return;

    for ( idx_t i = 0; i < N; ++i ) {
        if ( src_fields[i].dirty() ) { source().haloExchange( const_cast<Field&>( src_fields[i] ) ); }
    }

    array::DataType datatype = src_fields[0].datatype();
    int rank                 = src_fields[0].rank();

    ASSERT( rank > 1 );

    for ( idx_t i = 0; i < N; ++i ) {
        ASSERT( src_fields[i].datatype() == datatype );
        ASSERT( src_fields[i].rank() == rank );
        ASSERT( tgt_fields[i].datatype() == datatype );
    }

    if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 2 ) {
        execute_impl<double, 2>( *kernel_, src_fields, tgt_fields );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 2 ) {
        execute_impl<float, 2>( *kernel_, src_fields, tgt_fields );
    }

    for ( idx_t i = 0; i < N; ++i ) {
        tgt_fields[i].set_dirty();
    }
}


template <typename Kernel>
template <typename Value, int Rank>
void StructuredInterpolation3D<Kernel>::execute_impl( const Kernel& kernel, const FieldSet& src_fields,
                                                      FieldSet& tgt_fields ) const {
    if ( functionspace::PointCloud tgt = target() ) {
        const idx_t N  = src_fields.size();
        idx_t out_npts = target_lonlat_.shape( 0 );

        auto ghost    = array::make_view<int, 1>( target_ghost_ );
        auto lonlat   = array::make_view<double, 2>( target_lonlat_ );
        auto vertical = array::make_view<double, 1>( target_vertical_ );

        std::vector<array::ArrayView<Value, Rank> > src_view;
        std::vector<array::ArrayView<Value, 1> > tgt_view;
        src_view.reserve( N );
        tgt_view.reserve( N );

        for ( idx_t i = 0; i < N; ++i ) {
            src_view.emplace_back( array::make_view<Value, Rank>( src_fields[i] ) );
            tgt_view.emplace_back( array::make_view<Value, 1>( tgt_fields[i] ) );
        }

        constexpr NormaliseLongitude normalise( 0., 360. );
        atlas_omp_parallel {
            typename Kernel::Stencil stencil;
            typename Kernel::Weights weights;
            atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                if ( not ghost( n ) ) {
                    double x = normalise( lonlat( n, LON ) );
                    double y = lonlat( n, LAT );
                    double z = vertical( n );
                    kernel.compute_stencil( x, y, z, stencil );
                    kernel.compute_weights( x, y, z, stencil, weights );
                    for ( idx_t i = 0; i < N; ++i ) {
                        kernel.interpolate( stencil, weights, src_view[i], tgt_view[i], n );
                    }
                }
            }
        }
    }
    else if ( target_3d_ ) {
        const idx_t N  = src_fields.size();
        idx_t out_npts = target_3d_.shape( 0 );
        idx_t out_nlev = target_3d_.shape( 1 );

        const auto coords = array::make_view<double, 3, array::Intent::ReadOnly>( target_3d_ );

        std::vector<array::ArrayView<Value, Rank, array::Intent::ReadOnly> > src_view;
        std::vector<array::ArrayView<Value, Rank> > tgt_view;
        src_view.reserve( N );
        tgt_view.reserve( N );

        for ( idx_t i = 0; i < N; ++i ) {
            src_view.emplace_back( array::make_view<Value, Rank, array::Intent::ReadOnly>( src_fields[i] ) );
            tgt_view.emplace_back( array::make_view<Value, Rank>( tgt_fields[i] ) );
        }

        constexpr NormaliseLongitude normalise( 0., 360. );
        atlas_omp_parallel {
            typename Kernel::Stencil stencil;
            typename Kernel::Weights weights;
            atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                for ( idx_t k = 0; k < out_nlev; ++k ) {
                    double x = normalise( coords( n, k, LON ) );
                    double y = coords( n, k, LAT );
                    double z = coords( n, k, ZZ );
                    kernel.compute_stencil( x, y, z, stencil );
                    kernel.compute_weights( x, y, z, stencil, weights );
                    for ( idx_t i = 0; i < N; ++i ) {
                        kernel.interpolate( stencil, weights, src_view[i], tgt_view[i], n, k );
                    }
                }
            }
        }
    }
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
