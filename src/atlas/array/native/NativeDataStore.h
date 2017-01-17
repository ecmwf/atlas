/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/internals/atlas_defines.h"
#include "atlas/array/ArrayUtil.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace native {

template<typename Value>
class DataStore : public ArrayDataStore {
public:

    DataStore(size_t size) :
        data_store_(size) {
    }

    void clone_to_device() const {
    }

    void clone_from_device() const {
    }

    bool valid() const {
        return true;
    }

    void sync() const {
    }

    bool is_on_host() const {
        return true;
    }

    bool is_on_device() const {
        return false;
    }

    void reactivate_device_write_views() const {
    }

    void reactivate_host_write_views() const {
    }

    void* void_data_store() {
        return static_cast<void*>( &data_store_.front() );
    }

private:
    std::vector<Value> data_store_;
};

//------------------------------------------------------------------------------

template<typename Value>
class WrappedDataStore : public ArrayDataStore {
public:

    WrappedDataStore(Value* data_store) :
        data_store_(data_store) {
    }

    void clone_to_device() const {
    }

    void clone_from_device() const {
    }

    bool valid() const {
        return true;
    }

    void sync() const {
    }

    bool is_on_host() const {
        return true;
    }

    bool is_on_device() const {
        return false;
    }

    void reactivate_device_write_views() const {
    }

    void reactivate_host_write_views() const {
    }

    void* void_data_store() {
        return static_cast<void*>( data_store_ );
    }

private:
    Value* data_store_;
};

} // namespace native
} // namespace array
} // namespace atlas

