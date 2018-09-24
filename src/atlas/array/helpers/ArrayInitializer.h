/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>

#include "eckit/exception/Exceptions.h"

#include "atlas/array.h"
#include "atlas/array/DataType.h"
#include "atlas/array_fwd.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace helpers {

//------------------------------------------------------------------------------

template <unsigned int Rank>
struct array_initializer;

template <unsigned int PartDim>
struct array_initializer_partitioned;

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, unsigned int Dim>
struct array_initializer_impl {
    static void apply( Array const& orig, Array& array_resized ) {
        array_initializer_impl<Value, Rank, Dim>::apply( make_view<Value, Rank>( orig ),
                                                         make_view<Value, Rank>( array_resized ) );
    }

    template <typename... DimIndex>
    static void apply( ArrayView<Value, Rank> const&& orig, ArrayView<Value, Rank>&& array_resized, DimIndex... idxs ) {
        const size_t N = std::min( array_resized.shape( Dim ), orig.shape( Dim ) );
        for ( size_t i = 0; i < N; ++i ) {
            array_initializer_impl<Value, Rank, Dim + 1>::apply( std::move( orig ), std::move( array_resized ), idxs...,
                                                                 i );
        }
    }
};

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank>
struct array_initializer_impl<Value, Rank, Rank> {
    template <typename... DimIndex>
    static void apply( ArrayView<Value, Rank> const&& orig, ArrayView<Value, Rank>&& array_resized, DimIndex... idxs ) {
        array_resized( idxs... ) = orig( idxs... );
    }
};

//------------------------------------------------------------------------------

template <unsigned int Rank>
struct array_initializer {
    static void apply( Array const& orig, Array& array_resized ) {
        switch ( orig.datatype().kind() ) {
            case DataType::KIND_REAL64:
                return array_initializer_impl<double, Rank, 0>::apply( orig, array_resized );
            case DataType::KIND_REAL32:
                return array_initializer_impl<float, Rank, 0>::apply( orig, array_resized );
            case DataType::KIND_INT32:
                return array_initializer_impl<int, Rank, 0>::apply( orig, array_resized );
            case DataType::KIND_INT64:
                return array_initializer_impl<long, Rank, 0>::apply( orig, array_resized );
            case DataType::KIND_UINT64:
                return array_initializer_impl<unsigned long, Rank, 0>::apply( orig, array_resized );
            default: {
                std::stringstream err;
                err << "data kind " << orig.datatype().kind() << " not recognised.";
                throw eckit::BadParameter( err.str(), Here() );
            }
        }
    }
};

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, unsigned int Dim, unsigned int PartDim>
struct array_initializer_partitioned_val_impl {
    static void apply( Array const& orig, Array& dest, unsigned int pos, unsigned int offset ) {
        array_initializer_partitioned_val_impl<Value, Rank, Dim, PartDim>::apply(
            make_view<Value, Rank>( orig ), make_view<Value, Rank>( dest ), pos, offset );
    }

    template <typename... DimIndexPair>
    static void apply( ArrayView<Value, Rank> const&& orig, ArrayView<Value, Rank>&& dest, unsigned int pos,
                       unsigned int offset, DimIndexPair... idxs ) {
        for ( size_t i = 0; i < orig.shape( Dim ); ++i ) {
            unsigned int displ = i;
            if ( Dim == PartDim && i >= pos ) { displ += offset; }
            std::pair<int, int> pair_idx{i, displ};
            array_initializer_partitioned_val_impl<Value, Rank, Dim + 1, PartDim>::apply(
                std::move( orig ), std::move( dest ), pos, offset, idxs..., pair_idx );
        }
    }
};

// template< typename stdarray >
// inline std::string print_array(const stdarray& v)
// {
//   std::stringstream s;
//   s << "[ ";
//   for( int j=0; j<v.size(); ++ j ) {
//     s << v[j];
//     if( j != v.size()-1 ) s << " , ";
//   }
//   s << " ]";
//   return s.str();
// }

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, unsigned int PartDim>
struct array_initializer_partitioned_val_impl<Value, Rank, Rank, PartDim> {
    template <typename... DimIndexPair>
    static void apply( ArrayView<Value, Rank> const&& orig, ArrayView<Value, Rank>&& dest, unsigned int pos,
                       unsigned int offset, DimIndexPair... idxs ) {
        // Log::info() << print_array(std::array<int,Rank>{std::get<0>(idxs)...}) <<
        // " --> " << print_array(std::array<int,Rank>{std::get<1>(idxs)...}) << "
        // " <<  orig(std::get<0>(idxs)...) << std::endl;
        dest( std::get<1>( idxs )... ) = orig( std::get<0>( idxs )... );
    }
};

//------------------------------------------------------------------------------

template <unsigned int Rank, unsigned int PartDim>
struct array_initializer_partitioned_impl {
    static void apply( Array const& orig, Array& dest, unsigned int pos, unsigned int offset ) {
        switch ( orig.datatype().kind() ) {
            case DataType::KIND_REAL64:
                return array_initializer_partitioned_val_impl<double, Rank, 0, PartDim>::apply( orig, dest, pos,
                                                                                                offset );
            case DataType::KIND_REAL32:
                return array_initializer_partitioned_val_impl<float, Rank, 0, PartDim>::apply( orig, dest, pos,
                                                                                               offset );
            case DataType::KIND_INT32:
                return array_initializer_partitioned_val_impl<int, Rank, 0, PartDim>::apply( orig, dest, pos, offset );
            case DataType::KIND_INT64:
                return array_initializer_partitioned_val_impl<long, Rank, 0, PartDim>::apply( orig, dest, pos, offset );
            case DataType::KIND_UINT64:
                return array_initializer_partitioned_val_impl<unsigned long, Rank, 0, PartDim>::apply( orig, dest, pos,
                                                                                                       offset );
            default: {
                std::stringstream err;
                err << "data kind " << orig.datatype().kind() << " not recognised.";
                throw eckit::BadParameter( err.str(), Here() );
            }
        }
    }
};

//------------------------------------------------------------------------------

template <unsigned int PartDim>
struct array_initializer_partitioned {
    static void apply( Array const& orig, Array& dest, unsigned int pos, unsigned int offset ) {
        switch ( orig.rank() ) {
            case 1:
                return array_initializer_partitioned_impl<1, PartDim>::apply( orig, dest, pos, offset );
            case 2:
                return array_initializer_partitioned_impl<2, PartDim>::apply( orig, dest, pos, offset );
            case 3:
                return array_initializer_partitioned_impl<3, PartDim>::apply( orig, dest, pos, offset );
            case 4:
                return array_initializer_partitioned_impl<4, PartDim>::apply( orig, dest, pos, offset );
            case 5:
                return array_initializer_partitioned_impl<5, PartDim>::apply( orig, dest, pos, offset );
            case 6:
                return array_initializer_partitioned_impl<6, PartDim>::apply( orig, dest, pos, offset );
            case 7:
                return array_initializer_partitioned_impl<7, PartDim>::apply( orig, dest, pos, offset );
            case 8:
                return array_initializer_partitioned_impl<8, PartDim>::apply( orig, dest, pos, offset );
            case 9:
                return array_initializer_partitioned_impl<9, PartDim>::apply( orig, dest, pos, offset );
            default: {
                std::stringstream err;
                err << "too high Rank";
                throw eckit::BadParameter( err.str(), Here() );
            }
        }
    }
};

//------------------------------------------------------------------------------

}  // namespace helpers
}  // namespace array
}  // namespace atlas
