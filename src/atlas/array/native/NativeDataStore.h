/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/array/ArrayUtil.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace native {

template <typename Value>
class DataStore : public ArrayDataStore {
public:
    DataStore( size_t size ) : data_store_( size ) {}

    void cloneToDevice() const {}

    void cloneFromDevice() const {}

    bool valid() const { return true; }

    void syncHostDevice() const {}

    bool hostNeedsUpdate() const { return false; }

    bool deviceNeedsUpdate() const { return false; }

    void reactivateDeviceWriteViews() const {}

    void reactivateHostWriteViews() const {}

    void* voidDataStore() { return static_cast<void*>( &data_store_.front() ); }

    void* voidHostData() { return static_cast<void*>( &data_store_.front() ); }

    void* voidDeviceData() { return static_cast<void*>( &data_store_.front() ); }

private:
    std::vector<Value> data_store_;
};

//------------------------------------------------------------------------------

template <typename Value>
class WrappedDataStore : public ArrayDataStore {
public:
    WrappedDataStore( Value* data_store ) : data_store_( data_store ) {}

    void cloneToDevice() const {}

    void cloneFromDevice() const {}

    bool valid() const { return true; }

    void syncHostDevice() const {}

    bool hostNeedsUpdate() const { return true; }

    bool deviceNeedsUpdate() const { return false; }

    void reactivateDeviceWriteViews() const {}

    void reactivateHostWriteViews() const {}

    void* voidDataStore() { return static_cast<void*>( data_store_ ); }

    void* voidHostData() { return static_cast<void*>( data_store_ ); }

    void* voidDeviceData() { return static_cast<void*>( data_store_ ); }

private:
    Value* data_store_;
};

}  // namespace native
}  // namespace array
}  // namespace atlas
