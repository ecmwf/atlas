/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



/// @file ArrayView_iterator.h
/// This file contains the ArrayView_iterator and ArrayView_const_iterator class

#ifndef atlas_ArrayView_h
#error This file has to be included by atlas_ArrayView_h only
#endif

#ifndef atlas_ArrayView_iterator_h
#define atlas_ArrayView_iterator_h

#include <iterator>

//------------------------------------------------------------------------------------------------------

namespace atlas {

template <typename DATA_TYPE, int RANK=0> class ArrayView_iterator;
template <typename DATA_TYPE, int RANK=0> class ArrayView_const_iterator;

template <typename DATA_TYPE, int RANK>
class ArrayView_iterator : public std::iterator<std::forward_iterator_tag, DATA_TYPE>
{
    friend class ArrayView_const_iterator<DATA_TYPE,RANK>;
private:
  DATA_TYPE* p_;                    //!< raw pointer
  ArrayView<DATA_TYPE,RANK>* arr_;  //!< array to be iterated
  ArrayIdx loc_;                    //!< current position in array
  int fastest_idx_;                 //!< store fastest-moving index

public:

  /// @brief constructor (to be used for end of range)
  ArrayView_iterator();

  /// @brief constructor (to be used for begin of range)
  ArrayView_iterator(ArrayView<DATA_TYPE,RANK>* arr);

  /// @brief constructor from array at specific point
  ArrayView_iterator(ArrayView<DATA_TYPE,RANK>* arr, const ArrayIdx& loc);

  /// @brief constructor from other iterator
  ArrayView_iterator(const ArrayView_iterator& it);

  /// @brief pre-increment operator
  ArrayView_iterator& operator++()    { return increment(fastest_idx_); }

  /// @brief post-increment operator
  ArrayView_iterator operator++(int)  { ArrayView_iterator<DATA_TYPE,RANK> tmp(*this); operator++(); return tmp; }

  /// @brief equals operator
  bool operator==(const ArrayView_iterator& rhs) { return p_==rhs.p_; }

  /// @brief not-equals operator
  bool operator!=(const ArrayView_iterator& rhs) { return p_!=rhs.p_; }

  /// @brief dereference operator
  DATA_TYPE& operator*()      { return *p_; }

  /// @brief current position in array
  const ArrayIdx& pos() const { return loc_; }

private:
  ArrayView_iterator& increment(int d);

};

template <typename DATA_TYPE, int RANK>
class ArrayView_const_iterator : public std::iterator<std::forward_iterator_tag, DATA_TYPE>
{
  friend class ArrayView_iterator<DATA_TYPE,RANK>;
private:
  const DATA_TYPE* p_;                    //!< raw pointer
  const ArrayView<DATA_TYPE,RANK>* arr_;  //!< array to be iterated
  ArrayIdx loc_;                    //!< current position in array
  int fastest_idx_;                 //!< store fastest-moving index

public:

  /// @brief constructor (to be used for end of range)
  ArrayView_const_iterator();

  /// @brief constructor (to be used for begin of range)
  ArrayView_const_iterator(const ArrayView<DATA_TYPE,RANK>* arr);

  /// @brief constructor from array at specific point
  ArrayView_const_iterator(const ArrayView<DATA_TYPE,RANK>* arr, const ArrayIdx& loc);

  /// @brief constructor from other iterator
  ArrayView_const_iterator(const ArrayView_iterator<DATA_TYPE,RANK>& it);

  /// @brief constructor from other iterator
  ArrayView_const_iterator(const ArrayView_const_iterator& it);

  /// @brief pre-increment operator
  ArrayView_const_iterator& operator++()    { return increment(fastest_idx_); }

  /// @brief post-increment operator
  ArrayView_const_iterator operator++(int)  { ArrayView_const_iterator<DATA_TYPE,RANK> tmp(*this); operator++(); return tmp; }

  /// @brief equals operator
  bool operator==(const ArrayView_iterator<DATA_TYPE,RANK>& rhs) { return p_==rhs.p_; }

  /// @brief equals operator
  bool operator==(const ArrayView_const_iterator& rhs) { return p_==rhs.p_; }

  /// @brief not-equals operator
  bool operator!=(const ArrayView_iterator<DATA_TYPE,RANK>& rhs) { return p_!=rhs.p_; }

  /// @brief not-equals operator
  bool operator!=(const ArrayView_const_iterator& rhs) { return p_!=rhs.p_; }

  /// @brief dereference operator
  const DATA_TYPE& operator*() const { return *p_; }

  /// @brief current position in array
  const ArrayIdx& pos() const { return loc_; }

private:
  ArrayView_const_iterator& increment(int d);
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
