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

#include "atlas/array/ArrayUtil.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace gridtools {

template<typename gt_DataStore>
struct GridToolsDataStore : ArrayDataStore
{
    explicit GridToolsDataStore(gt_DataStore const *ds) : data_store_(ds) {}

    ~GridToolsDataStore() {
        assert(data_store_);
        delete data_store_;
    }

    void cloneToDevice() const {
        data_store_->clone_to_device();
    }

    void cloneFromDevice() const {
        data_store_->clone_from_device();
    }

    bool valid() const {
        return data_store_->valid();
    }

    void syncHostDevice() const {
        data_store_->sync();
    }

    bool hostNeedsUpdate() const {
        return data_store_->host_needs_update();
    }

    bool deviceNeedsUpdate() const {
        return data_store_->device_needs_update();
    }

    void reactivateDeviceWriteViews() const {
        data_store_->reactivate_device_write_views();
    }

    void reactivateHostWriteViews() const {
        data_store_->reactivate_host_write_views();
    }

    void* voidDataStore() {
        return static_cast<void*>(const_cast<gt_DataStore*>(data_store_));
    }

private:
    gt_DataStore const* data_store_;
};

} // namespace gridtools
} // namespace array
} // namespace atlas
