/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an size_tergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef atlas_ArrayUtil_h
#define atlas_ArrayUtil_h

#include <stddef.h>
#include <vector>

//------------------------------------------------------------------------------------------------------

namespace atlas {

template<typename T> struct remove_const          { typedef T type; };
template<typename T> struct remove_const<T const> { typedef T type; };

template<typename T> struct add_const          { typedef const typename remove_const<T>::type type; };
template<typename T> struct add_const<T const> { typedef const T type; };

typedef std::vector<size_t> ArrayShape;
typedef std::vector<size_t> ArrayStrides;
typedef std::vector<size_t> ArrayIdx;

inline ArrayShape make_shape() { return std::vector<size_t>(); }
inline ArrayShape make_shape(size_t size1) { return std::vector<size_t>(1,size1); }
inline ArrayShape make_shape(size_t size1, size_t size2) { std::vector<size_t> v(2); v[0]=size1; v[1]=size2; return v; }
inline ArrayShape make_shape(size_t size1, size_t size2, size_t size3) { std::vector<size_t> v(3); v[0]=size1; v[1]=size2; v[2]=size3; return v; }
inline ArrayShape make_shape(size_t size1, size_t size2, size_t size3, size_t size4) { std::vector<size_t> v(4); v[0]=size1; v[1]=size2; v[2]=size3; v[3]=size4; return v; }

inline ArrayStrides make_strides(size_t size1) { return std::vector<size_t>(1,size1); }
inline ArrayStrides make_strides(size_t size1, size_t size2) { std::vector<size_t> v(2); v[0]=size1; v[1]=size2; return v; }
inline ArrayStrides make_strides(size_t size1, size_t size2, size_t size3) { std::vector<size_t> v(3); v[0]=size1; v[1]=size2; v[2]=size3; return v; }
inline ArrayStrides make_strides(size_t size1, size_t size2, size_t size3, size_t size4) { std::vector<size_t> v(4); v[0]=size1; v[1]=size2; v[2]=size3; v[3]=size4; return v; }

inline ArrayIdx make_idx(size_t size1) { return std::vector<size_t>(1,size1); }
inline ArrayIdx make_idx(size_t size1, size_t size2) { std::vector<size_t> v(2); v[0]=size1; v[1]=size2; return v; }
inline ArrayIdx make_idx(size_t size1, size_t size2, size_t size3) { std::vector<size_t> v(3); v[0]=size1; v[1]=size2; v[2]=size3; return v; }
inline ArrayIdx make_idx(size_t size1, size_t size2, size_t size3, size_t size4) { std::vector<size_t> v(4); v[0]=size1; v[1]=size2; v[2]=size3; v[3]=size4; return v; }

class ArraySpec {
private:
  size_t size_;
  size_t rank_;
  ArrayShape shape_;
  ArrayStrides strides_;
  mutable std::vector<int> shapef_;
public:
  ArraySpec() : size_(0), rank_(0) {}
  ArraySpec( const ArrayShape& );
  size_t size() const { return size_; }
  size_t rank() const { return rank_; }
  const ArrayShape& shape() const { return shape_; }
  const ArrayStrides& strides() const { return strides_; }
  const std::vector<int>& shapef() const;
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
