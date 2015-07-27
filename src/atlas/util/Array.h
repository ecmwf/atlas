/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef atlas_Array_h
#define atlas_Array_h

#include <vector>

#include "eckit/memory/Owned.h"

#include "atlas/util/ArrayUtil.h"
#include "atlas/util/DataType.h"
#include "atlas/util/ArrayView.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

class ArrayBase : public eckit::Owned {
public:

  virtual std::string datatype() const = 0;

  void resize(const ArrayShape&);

  void resize(size_t size1);

  void resize(size_t size1, size_t size2);

  void resize(size_t size1, size_t size2, size_t size3);

  void resize(size_t size1, size_t size2, size_t size3, size_t size4);

  size_t size() const { return spec_.size(); }

  size_t rank() const { return spec_.rank(); }

  size_t stride(size_t i) const { return spec_.strides()[i]; }

  size_t shape(size_t i) const { return spec_.shape()[i]; }

  const ArrayStrides& strides() const { return spec_.strides(); }

  const ArrayShape& shape() const { return spec_.shape(); }

  const std::vector<int>& shapef() const { return spec_.shapef(); }

  /// @brief Access to raw data
  template <typename DATATYPE>       DATATYPE* data();
  template <typename DATATYPE> const DATATYPE* data() const;

private: // methods

  virtual void resize_data( size_t size )=0;

private:

  ArraySpec spec_;
};


template< typename DATA_TYPE >
class Array : public ArrayBase  {
public:

  typedef typename remove_const<DATA_TYPE>::type  value_type;
  typedef typename add_const<DATA_TYPE>::type     const_value_type;

public:

  Array() {}

  Array(const ArrayShape& shape) { resize(shape); }

  Array(size_t size)                                            { resize( make_shape(size) ); }

  Array(size_t size1, size_t size2)                             { resize( make_shape(size1,size2) ); }

  Array(size_t size1, size_t size2, size_t size3)               { resize( make_shape(size1,size2,size3) ); }

  Array(size_t size1, size_t size2, size_t size3, size_t size4) { resize( make_shape(size1,size2,size3,size4) ); }

public:

  const DATA_TYPE& operator()(size_t i) const                                         { return view_(i); }
        DATA_TYPE& operator()(size_t i)                                               { return view_(i); }
  const DATA_TYPE& operator()(size_t i, size_t j) const                               { return view_(i,j); }
        DATA_TYPE& operator()(size_t i, size_t j)                                     { return view_(i,j); }
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k) const                     { return view_(i,j,k); }
        DATA_TYPE& operator()(size_t i, size_t j, size_t k)                           { return view_(i,j,k); }
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l) const           { return view_(i,j,k,l); }
        DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l)                 { return view_(i,j,k,l); }
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const { return view_(i,j,k,l,m); }
        DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l, size_t m)       { return view_(i,j,k,l,m); }
  const DATA_TYPE& operator()(const ArrayIdx& idx) const                              { return view_(idx); }
        DATA_TYPE& operator()(const ArrayIdx& idx)                                    { return view_(idx); }

  virtual std::string datatype() const { return DataType::datatype<DATA_TYPE>(); }

  const DATA_TYPE& operator[](size_t i) const { return *(data()+i); }
        DATA_TYPE& operator[](size_t i)       { return *(data()+i); }

  const DATA_TYPE* data() const { return data_.data(); }
        DATA_TYPE* data()       { return data_.data(); }

  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<size(); ++n) data_[n]=scalar; }

  template< class InputIt >
  void assign( InputIt first, InputIt last ) { data_.assign(first,last); }

private:

  virtual void resize_data( size_t size );

private:

  std::vector<DATA_TYPE> data_;
  ArrayView<DATA_TYPE> view_;
};


template< typename DATA_TYPE>
void Array<DATA_TYPE>::resize_data( size_t size )
{
  data_.resize( size );
  view_ = ArrayView<DATA_TYPE>( *this );
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
