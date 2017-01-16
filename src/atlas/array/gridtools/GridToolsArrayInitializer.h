/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <vector>
#include "eckit/exception/Exceptions.h"
#include "atlas/array/DataType.h"
#include "atlas/array/array_fwd.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace gridtools {

//------------------------------------------------------------------------------

template <unsigned int RANK>
struct array_initializer;

template<unsigned int PartDim>
struct array_initializer_partitioned;

//------------------------------------------------------------------------------

template <typename Value, unsigned int RANK, unsigned int Dim>
struct array_initializer_impl {


  static void apply(Array const& orig, Array& array_resized) {
    array_initializer_impl<Value, RANK, Dim>::apply(
          make_view<Value, RANK>(orig),
          make_view<Value, RANK>(array_resized) );
  }


  template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& array_resized, DimIndex... idxs) {

    for(size_t i=0; i < orig.shape(Dim); ++i) {
      array_initializer_impl<Value, RANK, Dim+1>::apply(
            std::move(orig),
            std::move(array_resized),
            idxs...,
            i );
    }

  }


};

//------------------------------------------------------------------------------

template <typename Value, unsigned int RANK>
struct array_initializer_impl<Value, RANK, RANK> {

  template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig,ArrayView<Value, RANK>&& array_resized, DimIndex... idxs) {
      array_resized(idxs...) = orig(idxs...);
  }

};

//------------------------------------------------------------------------------

template <unsigned int RANK>
struct array_initializer {
  template <typename... DimIndex>
  static void apply(Array const& orig, Array& array_resized, DimIndex... idxs) {
    switch (orig.datatype().kind()) {
      case DataType::KIND_REAL64: return array_initializer_impl<double,        RANK, 0>::apply(orig, array_resized, idxs...);
      case DataType::KIND_REAL32: return array_initializer_impl<float,         RANK, 0>::apply(orig, array_resized, idxs...);
      case DataType::KIND_INT32:  return array_initializer_impl<int,           RANK, 0>::apply(orig, array_resized, idxs...);
      case DataType::KIND_INT64:  return array_initializer_impl<long,          RANK, 0>::apply(orig, array_resized, idxs...);
      case DataType::KIND_UINT64: return array_initializer_impl<unsigned long, RANK, 0>::apply(orig, array_resized, idxs...);
      default: {
        std::stringstream err;
        err << "data kind " << orig.datatype().kind() << " not recognised.";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }
};

//------------------------------------------------------------------------------

template<typename Value, unsigned int RANK, unsigned int Dim, unsigned int PartDim>
struct array_initializer_partitioned_val_impl {
  static void apply(Array const& orig, Array& dest, unsigned int pos, unsigned int offset) {
      auto view = make_view<Value, RANK>(orig);
      array_initializer_partitioned_val_impl<Value, RANK, Dim, PartDim>::apply(make_view<Value, RANK>(orig), make_view<Value, RANK>(dest), pos, offset);
  }
  template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& dest, unsigned int pos, unsigned int offset, DimIndex... idxs) {
      for(size_t i=0; i < orig.shape(Dim); ++i)
      {
          unsigned int displ = i;
          if(Dim == PartDim && i >= pos) {
              displ += offset;
          }
          std::pair<int,int> pair_idx{i,displ};
          array_initializer_partitioned_val_impl<Value, RANK, Dim+1, PartDim>::apply(std::move(orig), std::move(dest), pos, offset, idxs..., pair_idx);
      }
  }

};

//------------------------------------------------------------------------------

template <typename Value, unsigned int RANK, unsigned int PartDim>
struct array_initializer_partitioned_val_impl<Value, RANK, RANK, PartDim> {
    template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& dest, unsigned int pos, unsigned int offset, DimIndex... idxs) {
      dest(std::get<1>(idxs)...) = orig(std::get<0>(idxs)...);
  }
};

//------------------------------------------------------------------------------

template<unsigned int RANK, unsigned int PartDim>
struct array_initializer_partitioned_impl {
  static void apply( Array const& orig, Array& dest, unsigned int pos, unsigned int offset) {
      switch (orig.datatype().kind()) {
        case DataType::KIND_REAL64: return array_initializer_partitioned_val_impl<double, RANK, 0, PartDim>::apply(orig, dest, pos, offset);
        case DataType::KIND_REAL32: return array_initializer_partitioned_val_impl<float, RANK, 0, PartDim>::apply(orig, dest,pos, offset );
        case DataType::KIND_INT32:  return array_initializer_partitioned_val_impl<int, RANK, 0, PartDim>::apply(orig, dest, pos, offset);
        case DataType::KIND_INT64:  return array_initializer_partitioned_val_impl<long, RANK, 0, PartDim>::apply(orig, dest, pos, offset);
        case DataType::KIND_UINT64: return array_initializer_partitioned_val_impl<unsigned long, RANK, 0, PartDim>::apply(orig, dest, pos, offset);
        default: {
          std::stringstream err;
          err << "data kind " << orig.datatype().kind() << " not recognised.";
          throw eckit::BadParameter(err.str(), Here());
        }
      }
  }
};

//------------------------------------------------------------------------------

template<unsigned int PartDim>
struct array_initializer_partitioned {
  static void apply(Array const& orig, Array& dest, unsigned int pos, unsigned int offset) {
    switch (orig.rank()) {
      case 1: return array_initializer_partitioned_impl<1, PartDim>::apply(orig, dest, pos, offset);
      case 2: return array_initializer_partitioned_impl<2, PartDim>::apply(orig, dest, pos, offset);
      case 3: return array_initializer_partitioned_impl<3, PartDim>::apply(orig, dest, pos, offset);
      case 4: return array_initializer_partitioned_impl<4, PartDim>::apply(orig, dest, pos, offset);
      case 5: return array_initializer_partitioned_impl<5, PartDim>::apply(orig, dest, pos, offset);
      case 6: return array_initializer_partitioned_impl<6, PartDim>::apply(orig, dest, pos, offset);
      case 7: return array_initializer_partitioned_impl<7, PartDim>::apply(orig, dest, pos, offset);
      case 8: return array_initializer_partitioned_impl<8, PartDim>::apply(orig, dest, pos, offset);
      case 9: return array_initializer_partitioned_impl<9, PartDim>::apply(orig, dest, pos, offset);
      default: {
        std::stringstream err;
        err << "too high rank";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }
};

//------------------------------------------------------------------------------

} // namespace gridtools
} // namespace array
} // namespace atlas
