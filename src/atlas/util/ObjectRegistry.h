/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_ObjectRegistry_h
#define atlas_util_ObjectRegistry_h

#include <vector>
#include "eckit/exception/Exceptions.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace util {

template< typename T >
class ObjectRegistry
{
public:
  typedef size_t Id;
  Id add(const T&);
  T* get(Id) const;
  void remove(Id);
  static Id invalid() { return 0; }
private:
  std::vector<T*> store_;
};

//------------------------------------------------------------------------------------------------------

template<typename T>
typename ObjectRegistry<T>::Id ObjectRegistry<T>::add(const T& object)
{
  store_.push_back( const_cast<T*>(&object) );
  return store_.size();   // Id starts at 1
}

template<typename T>
T* ObjectRegistry<T>::get(Id id) const
{
  if( id == invalid() || id > store_.size() )
    throw eckit::Exception("Invalid id",Here());

  if( store_[id-1] == NULL )
    throw eckit::Exception("Trying to access object from registry, but it was previously removed");

  return store_[id-1];
}

template<typename T>
void ObjectRegistry<T>::remove(Id id)
{
  if( id == invalid() || id > store_.size() )
    throw eckit::Exception("Invalid id",Here());

  if( store_[id-1] == NULL )
    throw eckit::Exception("Trying to remove object from registry, but it was already removed");

  store_[id-1] = NULL;
}

//------------------------------------------------------------------------------------------------------

} // namespace util
} // namespace atlas

#endif
