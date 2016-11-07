#pragma once

#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

template <typename Value, unsigned int RANK, unsigned int Dim>
struct array_initializer_impl {

  static void apply(Array const& orig, Array& array_resized) {
      array_initializer_impl<Value, RANK, Dim>::apply(make_view<Value, RANK>(&orig), make_view<Value, RANK>(&array_resized));
  }
  template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& array_resized, DimIndex... idxs) {
      for(size_t i=0; i < orig.shape()[Dim]; ++i)
      {
          array_initializer_impl<Value, RANK, Dim+1>::apply(std::move(orig), std::move(array_resized), idxs..., i);
      }
  }
};

template <typename Value, unsigned int RANK>
struct array_initializer_impl<Value, RANK, RANK> {
    template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& array_resized, DimIndex... idxs) {
      array_resized(idxs...) = orig(idxs...);
  }
};

template <unsigned int RANK>
struct array_initializer {
  template <typename... DimIndex>
  static void apply(Array const& orig, Array& array_resized, DimIndex... idxs) {
    switch (orig.datatype().kind()) {
      case DataType::KIND_REAL64:
        array_initializer_impl<double, RANK, 0>::apply(orig, array_resized, idxs...);
        break;
      case DataType::KIND_REAL32:
        array_initializer_impl<float, RANK, 0>::apply(orig, array_resized, idxs...);
        break;
      case DataType::KIND_INT32:
        array_initializer_impl<int, RANK, 0>::apply(orig, array_resized, idxs...);
        break;
      case DataType::KIND_INT64:
        array_initializer_impl<long, RANK, 0>::apply(orig, array_resized, idxs...);
        break;
        ;
      default: {
        std::stringstream err;
        err << "data kind " << orig.datatype().kind() << " not recognised.";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }
};

}
}
#endif
