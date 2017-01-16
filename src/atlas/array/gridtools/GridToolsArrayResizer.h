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

#include "eckit/exception/Exceptions.h"

#include "atlas/array/array_fwd.h"

#include "atlas/array/ArrayUtil.h"
#include "atlas/array/DataType.h"
#include "atlas/array/ArrayView.h"

#include "atlas/array/gridtools/GridToolsArrayHelpers.h"
#include "atlas/array/gridtools/GridToolsArrayInitializer.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace gridtools {

template <int RANK> using UintSequence = ::gridtools::make_gt_integer_sequence<unsigned int, RANK>;

class GridToolsArrayResizer {
public:
  GridToolsArrayResizer( ArrayBase& array ) : array_(array) {}

  template <typename... Coords, typename = ::gridtools::all_integers<Coords...> >
  void resize(Coords... c) {
    if(sizeof...(c) != array_.rank()){
      std::stringstream err; err << "trying to resize an array of rank " << array_.rank() <<
                                    " by dimensions with rank " <<
                                    sizeof...(c) << std::endl;
      throw eckit::BadParameter(err.str(),Here());
    }

    check_dimension_lengths(array_.shape(), c...);

    if(!array_.data_store_->is_on_host()) {
        array_.data_store_->clone_from_device();
    }

    ArrayBase* array_resized = ArrayBase::create(array_.datatype(), ArrayShape{(unsigned int)c...});

    array_initializer<sizeof...(c)>::apply( array_, *array_resized);
    array_.data_store_.swap(array_resized->data_store_);
    array_.spec_ = array_resized->spec();

    delete array_resized;
  }

  void resize(const ArrayShape& shape)
  {
    assert(shape.size() > 0);
    switch (shape.size()) {
      case 1: return apply(shape, UintSequence<1>());
      case 2: return apply(shape, UintSequence<2>());
      case 3: return apply(shape, UintSequence<3>());
      case 4: return apply(shape, UintSequence<4>());
      case 5: return apply(shape, UintSequence<5>());
      case 6: return apply(shape, UintSequence<6>());
      case 7: return apply(shape, UintSequence<7>());
      case 8: return apply(shape, UintSequence<8>());
      case 9: return apply(shape, UintSequence<9>());
      default: {
        std::stringstream err;
        err << "shape not recognized";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }

private:

  template<typename UInt, UInt ... Indices>
  void apply(const ArrayShape& shape, ::gridtools::gt_integer_sequence<UInt, Indices...> ) {
    return resize(shape[Indices]...);
  }

  ArrayBase& array_;
};

//------------------------------------------------------------------------------

} // namespace gridtools
} // namespace array
} // namespace atlas
