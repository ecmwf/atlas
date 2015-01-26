/*
 * (C) Copyright 1996-2014 ECMWF.
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
#include "ArrayUtil.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

template< typename DATA_TYPE >
class Array {
public:
  typedef typename remove_const<DATA_TYPE>::type  value_type;
  typedef typename add_const<DATA_TYPE>::type     const_value_type;
public:
  Array() {}
  Array(int size) { resize( make_shape(size) ); }
  Array(int size1, int size2) { resize( make_shape(size1,size2) ); }
  Array(int size1, int size2, int size3) { resize( make_shape(size1,size2,size3) ); }
  Array(int size1, int size2, int size3, int size4) { resize( make_shape(size1,size2,size3,size4) ); }
  Array(const ArrayShape& shape) { resize(shape); }

  void resize(const ArrayShape& shape)
  {
    shape_= shape;
    strides_.resize(shape_.size());
    strides_[shape_.size()-1] = 1;
    for( int n=shape_.size()-2; n>=0; --n )
    {
      strides_[n] = strides_[n+1]*shape_[n+1];
    }
    std::size_t size = strides_[0]*shape_[0];
    data_.resize(size);
    rank_ = shape_.size();
  }
  const DATA_TYPE& operator[](int i) const { return *(data()+i); }
  DATA_TYPE&       operator[](int i)       { return *(data()+i); }
  std::size_t size() const { return data_.size(); }
  DATA_TYPE*       data() { return data_.data(); }
  const DATA_TYPE* data() const { return data_.data(); }
  int rank() const { return rank_; }
  int stride(int i) const { return strides_[i]; }
  int shape(int i) const { return shape_[i]; }
  const ArrayStrides& strides() const { return strides_; }
  const ArrayShape& shape() const { return shape_; }
  void operator=(const DATA_TYPE& scalar) { for(int n=0; n<size(); ++n) data_[n]=scalar; }
  template< class InputIt >
  void assign( InputIt first, InputIt last ) { data_.assign(first,last); }
private:
  int rank_;
  ArrayShape shape_;
  ArrayStrides strides_;
  std::vector<DATA_TYPE> data_;
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
