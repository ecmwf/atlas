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

#include <sstream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/array/DataType.h"
#include "atlas/array_fwd.h"
#include "atlas/runtime/Exception.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace helpers {

//------------------------------------------------------------------------------

struct array_initializer;

template <idx_t PartDim>
struct array_initializer_partitioned;

//------------------------------------------------------------------------------

template <typename Value, idx_t Rank, idx_t Dim>
struct array_initializer_impl {
    static void apply(Array const& orig, Array& array_resized) {
        array_initializer_impl<Value, Rank, Dim>::apply(make_view<Value, Rank>(orig),
                                                        make_view<Value, Rank>(array_resized));
    }

    template <typename... DimIndex>
    static void apply(ArrayView<const Value, Rank> const&& orig, ArrayView<Value, Rank>&& array_resized,
                      DimIndex... idxs) {
        const idx_t N = std::min(array_resized.shape(Dim), orig.shape(Dim));
        for (idx_t i = 0; i < N; ++i) {
            array_initializer_impl<Value, Rank, Dim + 1>::apply(std::move(orig), std::move(array_resized), idxs..., i);
        }
    }
};

//------------------------------------------------------------------------------

template <typename Value, idx_t Rank>
struct array_initializer_impl<Value, Rank, Rank> {
    template <typename... DimIndex>
    static void apply(ArrayView<const Value, Rank> const&& orig, ArrayView<Value, Rank>&& array_resized,
                      DimIndex... idxs) {
        array_resized(idxs...) = orig(idxs...);
    }
};

//------------------------------------------------------------------------------

struct array_initializer {

    static void apply(Array const& from, Array& to) {
        ATLAS_ASSERT(from.rank() == to.rank());
        switch (from.rank()) {
            case 1:
                apply_rank<1>(from, to);
                break;
            case 2:
                apply_rank<2>(from, to);
                break;
            case 3:
                apply_rank<3>(from, to);
                break;
            case 4:
                apply_rank<4>(from, to);
                break;
            case 5:
                apply_rank<5>(from, to);
                break;
            case 6:
                apply_rank<6>(from, to);
                break;
            case 7:
                apply_rank<7>(from, to);
                break;
            case 8:
                apply_rank<8>(from, to);
                break;
            case 9:
                apply_rank<9>(from, to);
                break;
            default:
                ATLAS_NOTIMPLEMENTED;
        }
    }

    template <idx_t Rank>
    static void apply_rank(Array const& orig, Array& array_resized) {
        switch (orig.datatype().kind()) {
            case DataType::KIND_REAL64:
                return array_initializer_impl<double, Rank, 0>::apply(orig, array_resized);
            case DataType::KIND_REAL32:
                return array_initializer_impl<float, Rank, 0>::apply(orig, array_resized);
            case DataType::KIND_INT32:
                return array_initializer_impl<int, Rank, 0>::apply(orig, array_resized);
            case DataType::KIND_INT64:
                return array_initializer_impl<long, Rank, 0>::apply(orig, array_resized);
            case DataType::KIND_UINT32:
                return array_initializer_impl<unsigned int, Rank, 0>::apply(orig, array_resized);
            case DataType::KIND_UINT64:
                return array_initializer_impl<unsigned long, Rank, 0>::apply(orig, array_resized);
            default: {
                std::stringstream err;
                err << "data kind " << orig.datatype().kind() << " not recognised.";
                throw_NotImplemented(err.str(), Here());
            }
        }
    }
};

//------------------------------------------------------------------------------

template <typename Value, idx_t Rank, idx_t Dim, idx_t PartDim>
struct array_initializer_partitioned_val_impl {
    static void apply(Array const& orig, Array& dest, idx_t pos, idx_t offset) {
        array_initializer_partitioned_val_impl<Value, Rank, Dim, PartDim>::apply(
            make_view<const Value, Rank>(orig), make_view<Value, Rank>(dest), pos, offset);
    }

    template <typename... DimIndexPair>
    static void apply(ArrayView<const Value, Rank>&& orig, ArrayView<Value, Rank>&& dest, idx_t pos, idx_t offset,
                      DimIndexPair... idxs) {
        for (idx_t i = 0; i < orig.shape(Dim); ++i) {
            idx_t displ = i;
            if (Dim == PartDim && i >= pos) {
                displ += offset;
            }
            std::pair<idx_t, idx_t> pair_idx{i, displ};
            array_initializer_partitioned_val_impl<Value, Rank, Dim + 1, PartDim>::apply(
                std::move(orig), std::move(dest), pos, offset, idxs..., pair_idx);
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

template <typename Value, idx_t Rank, idx_t PartDim>
struct array_initializer_partitioned_val_impl<Value, Rank, Rank, PartDim> {
    template <typename... DimIndexPair>
    static void apply(ArrayView<const Value, Rank>&& orig, ArrayView<Value, Rank>&& dest, idx_t /*pos*/,
                      idx_t /*offset*/, DimIndexPair... idxs) {
        // Log::info() << print_array(std::array<int,Rank>{std::get<0>(idxs)...}) <<
        // " --> " << print_array(std::array<int,Rank>{std::get<1>(idxs)...}) << "
        // " <<  orig(std::get<0>(idxs)...) << std::endl;
        dest(std::get<1>(idxs)...) = orig(std::get<0>(idxs)...);
    }
};

//------------------------------------------------------------------------------

template <idx_t Rank, idx_t PartDim>
struct array_initializer_partitioned_impl {
    static void apply(Array const& orig, Array& dest, idx_t pos, idx_t offset) {
        switch (orig.datatype().kind()) {
            case DataType::KIND_REAL64:
                return array_initializer_partitioned_val_impl<double, Rank, 0, PartDim>::apply(orig, dest, pos, offset);
            case DataType::KIND_REAL32:
                return array_initializer_partitioned_val_impl<float, Rank, 0, PartDim>::apply(orig, dest, pos, offset);
            case DataType::KIND_INT32:
                return array_initializer_partitioned_val_impl<int, Rank, 0, PartDim>::apply(orig, dest, pos, offset);
            case DataType::KIND_INT64:
                return array_initializer_partitioned_val_impl<long, Rank, 0, PartDim>::apply(orig, dest, pos, offset);
            case DataType::KIND_UINT32:
                return array_initializer_partitioned_val_impl<unsigned int, Rank, 0, PartDim>::apply(orig, dest, pos,
                                                                                                      offset);
            case DataType::KIND_UINT64:
                return array_initializer_partitioned_val_impl<unsigned long, Rank, 0, PartDim>::apply(orig, dest, pos,
                                                                                                      offset);
            default: {
                std::stringstream err;
                err << "data kind " << orig.datatype().kind() << " not recognised.";
                throw_NotImplemented(err.str(), Here());
            }
        }
    }
};

//------------------------------------------------------------------------------

template <idx_t PartDim>
struct array_initializer_partitioned {
    static void apply(const Array& orig, Array& dest, idx_t pos, idx_t offset) {
        switch (orig.rank()) {
            case 1:
                return array_initializer_partitioned_impl<1, PartDim>::apply(orig, dest, pos, offset);
            case 2:
                return array_initializer_partitioned_impl<2, PartDim>::apply(orig, dest, pos, offset);
            case 3:
                return array_initializer_partitioned_impl<3, PartDim>::apply(orig, dest, pos, offset);
            case 4:
                return array_initializer_partitioned_impl<4, PartDim>::apply(orig, dest, pos, offset);
            case 5:
                return array_initializer_partitioned_impl<5, PartDim>::apply(orig, dest, pos, offset);
            case 6:
                return array_initializer_partitioned_impl<6, PartDim>::apply(orig, dest, pos, offset);
            case 7:
                return array_initializer_partitioned_impl<7, PartDim>::apply(orig, dest, pos, offset);
            case 8:
                return array_initializer_partitioned_impl<8, PartDim>::apply(orig, dest, pos, offset);
            case 9:
                return array_initializer_partitioned_impl<9, PartDim>::apply(orig, dest, pos, offset);
            default: {
                std::stringstream err;
                err << "too high Rank";
                throw_NotImplemented(err.str(), Here());
            }
        }
    }
};

//------------------------------------------------------------------------------

}  // namespace helpers
}  // namespace array
}  // namespace atlas
