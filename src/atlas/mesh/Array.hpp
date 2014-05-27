/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef atlas_Array_hpp
#define atlas_Array_hpp

#include <vector>

//------------------------------------------------------------------------------------------------------

namespace atlas {

inline std::vector<int> Extents(int size1) { return std::vector<int>(1,size1); }
inline std::vector<int> Extents(int size1, int size2) { std::vector<int> extents(2); extents[0]=size1; extents[1]=size2; return extents; }
inline std::vector<int> Extents(int size1, int size2, int size3) { std::vector<int> extents(3); extents[0]=size1; extents[1]=size2; extents[2]=size3; return extents; }
inline std::vector<int> Extents(int size1, int size2, int size3, int size4) { std::vector<int> extents(4); extents[0]=size1; extents[1]=size2; extents[2]=size3; extents[3]=size4; return extents; }

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class Array {
public:
  Array() {}
  Array(int size) { resize( Extents(size) ); }
  Array(int size1, int size2) { resize( Extents(size1,size2) ); }
  Array(int size1, int size2, int size3) { resize( Extents(size1,size2,size3) ); }
  Array(int size1, int size2, int size3, int size4) { resize( Extents(size1,size2,size3,size4) ); }
  Array(const std::vector<int>& extents) { resize(extents); }
  
  void resize(const std::vector<int>& extents)
  {
    extents_= extents;
    strides_.resize(extents_.size());
    strides_[extents_.size()-1] = 1;
    for( int n=extents_.size()-2; n>=0; --n )
    {
      strides_[n] = strides_[n+1]*extents_[n+1];
    }
    std::size_t size = strides_[0]*extents_[0];
    data_.resize(size);
  }
  const DATA_TYPE& operator[](int i) const { return *(data()+i); }
  DATA_TYPE&       operator[](int i)       { return *(data()+i); }
  std::size_t size() const { return data_.size(); }
  DATA_TYPE*       data() { return data_.data(); }
  const DATA_TYPE* data() const { return data_.data(); }
  int stride(int i) const { return strides_[i]; }
  int extent(int i) const { return extents_[i]; }
  const std::vector<int>& strides() const { return strides_; }
  const std::vector<int>& extents() const { return extents_; }
  void operator=(const DATA_TYPE& scalar) { for(int n=0; n<size(); ++n) data_[n]=scalar; }
private:
  std::vector<int> extents_;
  std::vector<int> strides_;
  std::vector<DATA_TYPE> data_;
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
