/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>

#include "atlas/array/ArrayIdx.h"
#include "atlas/array/ArrayLayout.h"
#include "atlas/array/ArrayShape.h"
#include "atlas/array/ArraySpec.h"
#include "atlas/array/ArrayStrides.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <typename T>
struct remove_const {
    typedef T type;
};
template <typename T>
struct remove_const<T const> {
    typedef T type;
};

template <typename T>
struct add_const {
    typedef const typename remove_const<T>::type type;
};
template <typename T>
struct add_const<T const> {
    typedef const T type;
};
#endif

class ArrayDataStore {
public:
    virtual ~ArrayDataStore() {}
    virtual void updateDevice() const               = 0;
    virtual void updateHost() const                 = 0;
    virtual bool valid() const                      = 0;
    virtual void syncHostDevice() const             = 0;
    virtual void allocateDevice() const             = 0;
    virtual void deallocateDevice() const           = 0;
    virtual bool deviceAllocated() const            = 0;
    virtual bool hostNeedsUpdate() const            = 0;
    virtual bool deviceNeedsUpdate() const          = 0;
    virtual void setHostNeedsUpdate(bool) const     = 0;
    virtual void setDeviceNeedsUpdate(bool) const   = 0;
    virtual void reactivateDeviceWriteViews() const = 0;
    virtual void reactivateHostWriteViews() const   = 0;
    virtual void* voidDataStore()                   = 0;
    virtual void* voidHostData()                    = 0;
    virtual void* voidDeviceData()                  = 0;
    virtual void accMap() const                     = 0;
    virtual void accUnmap() const                   = 0;
    virtual bool accMapped() const                  = 0;
    template <typename Value>
    Value* hostData() {
        return static_cast<Value*>(voidHostData());
    }
    template <typename Value>
    Value* deviceData() {
        return static_cast<Value*>(voidDeviceData());
    }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <int Dim>
static constexpr char array_dim() {
    return Dim == 0 ? 'i' : (Dim == 1 ? 'j' : (Dim == 2 ? 'k' : (Dim == 3 ? 'l' : (Dim == 4 ? 'm' : ('*')))));
}

void throw_OutOfRange(const std::string& class_name, char idx_str, int idx, int max);
#endif

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
