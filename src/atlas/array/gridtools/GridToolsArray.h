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

#include "atlas/array/ArrayUtil.h"
#include "atlas/array/DataType.h"
#include "atlas/array/Array.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

class Array;

template <typename Value>
class ArrayT;

template <typename Value>
class ArrayCreatorT
{
public:
  Array* create(size_t dim0);
  Array* create(size_t dim0, size_t dim1);
  Array* create(size_t dim0, size_t dim1, size_t dim2);
  Array* create(size_t dim0, size_t dim1, size_t dim2, size_t dim3);
  Array* create(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4);
  Array* create( const ArrayShape& );
  Array* create( const ArrayShape&, const ArrayLayout& );
  Array* wrap(Value* data, const ArrayShape& );
  Array* wrap(Value* data, const ArraySpec& );

public:
  Array* create(ArrayDataStore*, const ArraySpec&);
};


//------------------------------------------------------------------------------

class Array : public ArrayBase {
public:

  static Array* create( array::DataType, const ArrayShape& );

  static Array* create( array::DataType, const ArrayShape&, const ArrayLayout& );

  template<typename Value> static Array* create(size_t dim0) {
    return ArrayCreatorT<Value>::create(dim0);
  }

  template<typename Value> static Array* create(size_t dim0, size_t dim1) {
    return ArrayCreatorT<Value>::create(dim0,dim1);
  }

  template<typename Value> static Array* create(size_t dim0, size_t dim1, size_t dim2) {
    return ArrayCreatorT<Value>::create(dim0,dim1,dim2);
  }

  template<typename Value> static Array* create(size_t dim0, size_t dim1, size_t dim2, size_t dim3) {
    return ArrayCreatorT<Value>::create(dim0,dim1,dim2,dim3);
  }

  template<typename Value> static Array* create(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
    return ArrayCreatorT<Value>::create(dim0,dim1,dim2,dim3,dim4);
  }

  template<typename Value> static Array* create(const ArrayShape& shape) {
    return create( array::DataType::create<Value>(), shape );
  }

  template<typename Value> static Array* create(const ArrayShape& shape, const ArrayLayout& layout) {
    return create( array::DataType::create<Value>(), shape, layout );
  }

  template <typename Value> static Array* wrap(Value* data, const ArrayShape& shape) {
    return ArrayCreatorT<Value>::wrap(data,shape);
  }

  template <typename Value> static Array* wrap(Value* data, const ArraySpec& spec) {
    return ArrayCreatorT<Value>::wrap(data,spec);
  }

};

//------------------------------------------------------------------------------

template<typename DATA_TYPE>
class ArrayT : public Array  {

  template <typename T>
  friend class ArrayTConstructor;
  friend class GridToolsArrayResizer;
  friend class GridToolsArrayInsert;

public:

  // These constructors are used to create an ArrayT on the stack.
  ArrayT(size_t dim0);
  ArrayT(size_t dim0, size_t dim1);
  ArrayT(size_t dim0, size_t dim1, size_t dim2);
  ArrayT(size_t dim0, size_t dim1, size_t dim2, size_t dim3);
  ArrayT(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4);

  ArrayT(const ArrayShape&);
  ArrayT(const ArrayShape&, const ArrayLayout&);

  ArrayT(const ArraySpec&);

public:

  virtual void insert(size_t idx1, size_t size1);

  virtual void resize(const ArrayShape&);
  virtual void resize(size_t dim0);
  virtual void resize(size_t dim0, size_t dim1);
  virtual void resize(size_t dim0, size_t dim1, size_t dim2);
  virtual void resize(size_t dim0, size_t dim1, size_t dim2, size_t dim3);
  virtual void resize(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4);

  virtual array::DataType datatype() const { return array::DataType::create<DATA_TYPE>(); }

  virtual size_t sizeof_data() const {return sizeof(DATA_TYPE);}

  virtual void dump(std::ostream& os) const;

private:
  // This constructor is used through the Array::create() or the Array::wrap() methods
  template <typename T>
  friend class ArrayCreatorT;
  ArrayT(ArrayDataStore*, const ArraySpec&);
};

//------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
