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

//------------------------------------------------------------------------------

class ArrayBackendInsert {

public:
  ArrayBackendInsert( Array& array ) : array_(array) {}

  void insert(size_t idx1, size_t size1)
  {
      if(!array_.is_on_host()) {
          array_.clone_from_device();
      }

      ArrayShape nshape = array_.shape();
      if(idx1 > nshape[0]) {
          throw eckit::BadParameter("can not insert into an array at a position beyond its size", Here());
      }
      nshape[0] += size1;

      Array* array_resized = Array::create(array_.datatype(), nshape);

      array_initializer_partitioned<0>::apply( array_, *array_resized, idx1, size1);
      array_.data_store_.swap(array_resized->data_store_);
      array_.spec_ = array_resized->spec();

      delete array_resized;
  }

private:

  Array& array_;
};

//------------------------------------------------------------------------------

} // namespace gridtools
} // namespace array
} // namespace atlas
