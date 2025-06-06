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

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {

template <typename Kernel>
double StructuredInterpolation3D<Kernel>::convert_units_multiplier( const Field& field ) {
    std::string units = field.metadata().getString( "units", "degrees" );
    if ( units == "degrees" ) {
        return 1.;
    }
    if ( units == "radians" ) {
        return 180. / M_PI;
    }
    ATLAS_NOTIMPLEMENTED;
}

template <typename Kernel>
StructuredInterpolation3D<Kernel>::StructuredInterpolation3D( const Method::Config& config ) :
    Method( config ),
    matrix_free_{false},
    limiter_{false} {
    config.get( "matrix_free", matrix_free_ );
    config.get( "limiter", limiter_ );

    if ( not matrix_free_ ) {
        throw_NotImplemented( "Matrix-free StructuredInterpolation3D not implemented", Here() );
    }
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::do_setup( const Grid& source, const Grid& target, const Cache& ) {
    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }
    ATLAS_ASSERT( StructuredGrid( source ) );
    FunctionSpace source_fs = functionspace::StructuredColumns( source, option::halo( kernel_->stencil_halo() ) );
    FunctionSpace target_fs = functionspace::PointCloud( target );

    do_setup( source_fs, target_fs );
}

template <typename Kernel>
void StructuredInterpolation3D<Kernel>::do_setup( const FunctionSpace& source, const FunctionSpace& target, const Cache& cache) {
    ATLAS_TRACE( "StructuredInterpolation3D<" + Kernel::className() + ">::do_setup(FunctionSpace source, FunctionSpace target, const Cache)" );
    if (! matrix_free_ && interpolation::MatrixCache(cache)) {
        setMatrix(cache);
        source_ = source;
        target_ = target;
        ATLAS_ASSERT(matrix().rows() == target.size());
        ATLAS_ASSERT(matrix().cols() == source.size());
        return;
    }
    else {
        do_setup( source, target );
    }
}

template <typename Kernel>
void StructuredInterpolation3D<Kernel>::do_setup( const FunctionSpace& source, const FunctionSpace& target ) {
    ATLAS_TRACE( "StructuredInterpolation<" + Kernel::className() + ">::do_setup()" );

    source_ = source;
    target_ = target;

    if ( functionspace::PointCloud tgt = target ) {
        target_lonlat_   = tgt.lonlat();
        target_vertical_ = tgt.vertical();
        target_ghost_    = tgt.ghost();
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }

    setup( source );
}

template <typename Kernel>
void StructuredInterpolation3D<Kernel>::do_setup( const FunctionSpace& source, const Field& target ) {
    ATLAS_TRACE( "StructuredInterpolation<" + Kernel::className() + ">::do_setup(FunctionSpace source, Field target)" );

    source_ = source;

    if ( target.functionspace() ) {
        target_ = target.functionspace();
    }

    target_3d_ = target;

    setup( source );
}

template <typename Kernel>
void StructuredInterpolation3D<Kernel>::do_setup( const FunctionSpace& source, const FieldSet& target ) {
    ATLAS_TRACE( "StructuredInterpolation<" + Kernel::className() + ">::do_setup(FunctionSpace source,FieldSet target)" );

    source_ = source;

    ATLAS_ASSERT( target.size() >= 3 );
    if ( target[0].functionspace() ) {
        target_ = target[0].functionspace();
    }

    target_xyz_ = target;

    setup( source );
}

template <typename Kernel>
void StructuredInterpolation3D<Kernel>::print( std::ostream& ) const {
    ATLAS_NOTIMPLEMENTED;
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::setup( const FunctionSpace& source ) {
    kernel_.reset( new Kernel( source, util::Config( "limiter", limiter_ ) ) );
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::do_execute( const Field& src_field, Field& tgt_field, Metadata& metadata ) const {
    FieldSet tgt( tgt_field );
    do_execute( FieldSet( src_field ), tgt, metadata );
}


template <typename Kernel>
void StructuredInterpolation3D<Kernel>::do_execute( const FieldSet& src_fields, FieldSet& tgt_fields, Metadata& metadata ) const {
    if ( not matrix_free_ ) {
        Method::do_execute( src_fields, tgt_fields, metadata );
        return;
    }

    const idx_t N = src_fields.size();
    ATLAS_ASSERT( N == tgt_fields.size() );

    if ( N == 0 )
        return;

    haloExchange( src_fields );

    ATLAS_TRACE( "StructuredInterpolation<" + Kernel::className() + ">::do_execute()" );
    
    array::DataType datatype = src_fields[0].datatype();
    int rank                 = src_fields[0].rank();

    ATLAS_ASSERT( rank > 1 );

    for ( idx_t i = 0; i < N; ++i ) {
        ATLAS_ASSERT( src_fields[i].datatype() == datatype );
        ATLAS_ASSERT( src_fields[i].rank() == rank );
        ATLAS_ASSERT( tgt_fields[i].datatype() == datatype );
    }

    if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 2 ) {
        execute_impl<double, 2>( *kernel_, src_fields, tgt_fields );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 2 ) {
        execute_impl<float, 2>( *kernel_, src_fields, tgt_fields );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 3 ) {
        execute_impl<double, 3>( *kernel_, src_fields, tgt_fields );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 3 ) {
        execute_impl<float, 3>( *kernel_, src_fields, tgt_fields );
    }

    tgt_fields.set_dirty();
}


template <typename Kernel>
template <typename Value, int Rank>
void StructuredInterpolation3D<Kernel>::execute_impl( const Kernel& kernel, const FieldSet& src_fields,
                                                      FieldSet& tgt_fields ) const {
    const idx_t N = src_fields.size();

    auto make_src_view = [&]( const FieldSet& src_fields ) {
        std::vector<array::ArrayView<const Value, Rank> > src_view;
        src_view.reserve( N );
        for ( idx_t i = 0; i < N; ++i ) {
            src_view.emplace_back( array::make_view<const Value, Rank>( src_fields[i] ) );
        }
        return src_view;
    };

    // Assertions
    ATLAS_ASSERT( tgt_fields.size() == src_fields.size() );
    idx_t tgt_rank = -1;
    for ( auto& f : tgt_fields ) {
        if ( tgt_rank == -1 )
            tgt_rank = f.rank();
        if ( f.rank() != tgt_rank ) {
            throw_Exception( "target fields don't all have the same rank!", Here() );
        }
    }

    if ( functionspace::PointCloud( target() ) && tgt_rank == 1 ) {
        const idx_t out_npts = target_lonlat_.shape( 0 );

        const auto ghost    = array::make_view<int, 1>( target_ghost_ );
        const auto lonlat   = array::make_view<double, 2>( target_lonlat_ );
        const auto vertical = array::make_view<double, 1>( target_vertical_ );

        auto src_view = make_src_view( src_fields );

        constexpr int TargetRank = 1;
        std::vector<array::ArrayView<Value, TargetRank> > tgt_view;
        tgt_view.reserve( N );
        for ( idx_t i = 0; i < N; ++i ) {
            tgt_view.emplace_back( array::make_view<Value, TargetRank>( tgt_fields[i] ) );
        }

        const double convert_units = convert_units_multiplier( target_lonlat_ );
        atlas_omp_parallel {
            typename Kernel::Stencil stencil;
            typename Kernel::Weights weights;
            atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                if ( not ghost( n ) ) {
                    double x = lonlat( n, LON ) * convert_units;
                    double y = lonlat( n, LAT ) * convert_units;
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
    else if ( target_3d_ && tgt_rank == Rank ) {
        const idx_t out_npts = target_3d_.shape( 0 );
        const idx_t out_nlev = target_3d_.shape( 1 );

        const auto coords   = array::make_view<const double, 3>( target_3d_ );
        const auto src_view = make_src_view( src_fields );

        constexpr int TargetRank = Rank;
        std::vector<array::ArrayView<Value, TargetRank> > tgt_view;
        tgt_view.reserve( N );

        for ( idx_t i = 0; i < N; ++i ) {
            tgt_view.emplace_back( array::make_view<Value, TargetRank>( tgt_fields[i] ) );

            if ( Rank == 3 &&
                 ( src_fields[i].stride( Rank - 1 ) != 1 || tgt_fields[i].stride( TargetRank - 1 ) != 1 ) ) {
                throw_Exception(
                    "Something will go seriously wrong if we continue from here as "
                    "the implementation assumes stride=1 for fastest moving index (variables).",
                    Here() );
            }
        }

        const double convert_units = convert_units_multiplier( target_3d_ );

        atlas_omp_parallel {
            typename Kernel::Stencil stencil;
            typename Kernel::Weights weights;
            atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                for ( idx_t k = 0; k < out_nlev; ++k ) {
                    double x = coords( n, k, LON ) * convert_units;
                    double y = coords( n, k, LAT ) * convert_units;
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
    else if ( not target_xyz_.empty() && tgt_rank == Rank ) {
        const idx_t out_npts = target_xyz_[0].shape( 0 );
        const idx_t out_nlev = target_xyz_[0].shape( 1 );

        const auto xcoords  = array::make_view<double, 2>( target_xyz_[LON] );
        const auto ycoords  = array::make_view<double, 2>( target_xyz_[LAT] );
        const auto zcoords  = array::make_view<double, 2>( target_xyz_[ZZ] );
        const auto src_view = make_src_view( src_fields );

        constexpr int TargetRank = Rank;
        std::vector<array::ArrayView<Value, TargetRank> > tgt_view;
        tgt_view.reserve( N );

        for ( idx_t i = 0; i < N; ++i ) {
            tgt_view.emplace_back( array::make_view<Value, TargetRank>( tgt_fields[i] ) );

            if ( Rank == 3 &&
                 ( src_fields[i].stride( Rank - 1 ) != 1 || tgt_fields[i].stride( TargetRank - 1 ) != 1 ) ) {
                throw_Exception(
                    "Something will go seriously wrong if we continue from here as "
                    "the implementation assumes stride=1 for fastest moving index (variables).",
                    Here() );
            }
        }

        const double convert_units = convert_units_multiplier( target_xyz_[LON] );

        atlas_omp_parallel {
            typename Kernel::Stencil stencil;
            typename Kernel::Weights weights;
            atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                for ( idx_t k = 0; k < out_nlev; ++k ) {
                    double x = xcoords( n, k ) * convert_units;
                    double y = ycoords( n, k ) * convert_units;
                    double z = zcoords( n, k );
                    kernel.compute_stencil( x, y, z, stencil );
                    kernel.compute_weights( x, y, z, stencil, weights );
                    for ( idx_t i = 0; i < N; ++i ) {
                        kernel.interpolate( stencil, weights, src_view[i], tgt_view[i], n, k );
                    }
                }
            }
        }
    }

    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
