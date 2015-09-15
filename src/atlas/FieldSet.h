/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @author Pedro Maciel
/// @date Jan 2015

#ifndef atlas_FieldSet_H
#define atlas_FieldSet_H

#include <vector>
#include <iterator>


#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/Field.h"

namespace atlas {

class FieldSet;

template <typename T> class SharedPtrVector_iterator;
template <typename T> class SharedPtrVector_const_iterator;

template< typename T>
class SharedPtrVector_iterator : public std::iterator<std::forward_iterator_tag, T>
{
  friend class SharedPtrVector_const_iterator<T>;
private:
  typedef std::vector< eckit::SharedPtr<T> >  SharedPtrVector;

  SharedPtrVector* arr_;  //!< array to be iterated  
  T* p_; //!< raw pointer

public:

  /// @brief constructor (to be used for end of range)
  SharedPtrVector_iterator() : arr_(0), p_(0) {}

  /// @brief constructor (to be used for begin of range)
  SharedPtrVector_iterator(SharedPtrVector* arr) : arr_(arr) { if( arr_->size() ) p_=(*arr_)[0].get(); }

  /// @brief constructor from other iterator
  SharedPtrVector_iterator(const SharedPtrVector_iterator<T>& it) : arr_(it.arr_), p_(it.p_) {}

  /// @brief pre-increment operator
  SharedPtrVector_iterator& operator++()    { return increment(1); }

  /// @brief post-increment operator
  SharedPtrVector_iterator operator++(int)  { SharedPtrVector_iterator tmp(*this); operator++(); return tmp; }

  /// @brief equals operator
  bool operator==(const SharedPtrVector_iterator& rhs) const { return p_==rhs.p_; }

  /// @brief not-equals operator
  bool operator!=(const SharedPtrVector_iterator& rhs) const { return p_!=rhs.p_; }

  /// @brief dereference operator
  T& operator*() { return *p_; }

  /// @brief dereference operator
  T* operator->() { return p_; }

private:
  SharedPtrVector_iterator& increment(size_t d)
  {
    p_ += d;
    if( d==arr_->size() )
      p_ = 0;
    return *this;
  }

};

template< typename T>
class SharedPtrVector_const_iterator : public std::iterator<std::forward_iterator_tag, T>
{
private:
  typedef const std::vector< eckit::SharedPtr<T> >  SharedPtrVector;

  SharedPtrVector* arr_;  //!< array to be iterated  
  T* p_; //!< raw pointer

public:

  /// @brief constructor (to be used for end of range)
  SharedPtrVector_const_iterator() : arr_(0), p_(0) {}

  /// @brief constructor (to be used for begin of range)
  SharedPtrVector_const_iterator(SharedPtrVector* arr) : arr_(arr) { if( arr_->size() ) p_=(*arr_)[0].get(); }

  /// @brief constructor from other iterator
  SharedPtrVector_const_iterator(const SharedPtrVector_const_iterator& it) : arr_(it.arr_), p_(it.p_) {}

  /// @brief constructor from other iterator
  SharedPtrVector_const_iterator(const SharedPtrVector_iterator<T>& it) : arr_(it.arr_), p_(it.p_) {}

  /// @brief pre-increment operator
  SharedPtrVector_const_iterator& operator++()    { return increment(1); }

  /// @brief post-increment operator
  SharedPtrVector_const_iterator operator++(int)  { SharedPtrVector_const_iterator tmp(*this); operator++(); return tmp; }

  /// @brief equals operator
  bool operator==(const SharedPtrVector_const_iterator& rhs) const { return p_==rhs.p_; }

  /// @brief not-equals operator
  bool operator!=(const SharedPtrVector_const_iterator& rhs) const { return p_!=rhs.p_; }

  /// @brief dereference operator
  const T& operator*() const { return *p_; }

  /// @brief dereference operator
  const T* operator->() const { return p_; }

private:
  SharedPtrVector_const_iterator& increment(size_t d)
  {
    p_ += d;
    if( d==arr_->size() )
      p_ = 0;
    return *this;
  }

};


/**
 * @brief Represents a set of fields, where order is preserved (no ownership)
 */
class FieldSet : public eckit::Owned {

public: // types

  typedef eckit::SharedPtr< FieldSet > Ptr;
  typedef SharedPtrVector_iterator<Field> iterator;
  typedef SharedPtrVector_const_iterator<Field> const_iterator;
public: // methods

  /// Constructs an empty FieldSet
  FieldSet(const std::string& name = "untitled");

  size_t size() const { return  fields_.size(); }
  bool empty()  const { return !fields_.size(); }

  const std::string& name() const { return name_; }
        std::string& name()       { return name_; }

  const Field& operator[](const size_t& i) const { return field(i); }
        Field& operator[](const size_t& i)       { return field(i); }

  const Field& field(const size_t& i) const { ASSERT(i<size()); return *fields_[i]; }
        Field& field(const size_t& i)       { ASSERT(i<size()); return *fields_[i]; }

  std::vector< std::string > field_names() const;

  void add(const Field&);

  bool has_field(const std::string& name) const;

  Field& field(const std::string& name) const;

  iterator begin() { return iterator(&fields_); }
  iterator end()   { return iterator(); }
  const_iterator begin() const { return const_iterator(&fields_); }
  const_iterator end() const   { return const_iterator(); }
  const_iterator cbegin() const { return const_iterator(&fields_); }
  const_iterator cend() const { return const_iterator(); }

protected: // data

  std::vector< eckit::SharedPtr<Field> >  fields_;  ///< field handle storage
  std::string                             name_;    ///< internal name
  std::map< std::string, size_t >         index_;   ///< name-to-index map, to refer fields by name
};


// C wrapper interfaces to C++ routines
extern "C"
{
  FieldSet* atlas__FieldSet__new           (char* name);
  void      atlas__FieldSet__delete        (FieldSet* This);
  void      atlas__FieldSet__add_field     (FieldSet* This, Field* field);
  int       atlas__FieldSet__has_field     (FieldSet* This, char* name);
  int       atlas__FieldSet__size          (FieldSet* This);
  Field*    atlas__FieldSet__field_by_name (FieldSet* This, char* name);
  Field*    atlas__FieldSet__field_by_idx  (FieldSet* This, int idx);
}


} // namespace atlas


#endif
