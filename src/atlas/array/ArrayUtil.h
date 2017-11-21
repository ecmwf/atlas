/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an size_tergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include "atlas/array/ArrayShape.h"
#include "atlas/array/ArrayStrides.h"
#include "atlas/array/ArrayLayout.h"
#include "atlas/array/ArrayIdx.h"
#include "atlas/array/ArraySpec.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

template<typename T> struct remove_const          { typedef T type; };
template<typename T> struct remove_const<T const> { typedef T type; };

template<typename T> struct add_const          { typedef const typename remove_const<T>::type type; };
template<typename T> struct add_const<T const> { typedef const T type; };

class ArrayDataStore
{
public:
  virtual ~ArrayDataStore() {}
  virtual void cloneToDevice() const = 0;
  virtual void cloneFromDevice() const = 0;
  virtual bool valid() const = 0;
  virtual void syncHostDevice() const = 0;
  virtual bool hostNeedsUpdate() const = 0;
  virtual bool deviceNeedsUpdate() const = 0;
  virtual void reactivateDeviceWriteViews() const = 0;
  virtual void reactivateHostWriteViews() const = 0;
  virtual void* voidDataStore() = 0;
  virtual void* voidHostData() = 0;
  virtual void* voidDeviceData() = 0;
  template <typename Value> Value* hostData()   { return static_cast<Value*>(voidHostData());   }
  template <typename Value> Value* deviceData() { return static_cast<Value*>(voidDeviceData()); }
};

template < int Dim >
static constexpr char array_dim() {
    return
        Dim == 0 ? 'i' :(
        Dim == 1 ? 'j' :(
        Dim == 2 ? 'k' :(
        Dim == 3 ? 'l' :(
        Dim == 4 ? 'm' :(
        '*')))));
}

void throw_OutOfRange( const std::string& class_name, char idx_str, int idx, int max );

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
