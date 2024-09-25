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

#include "atlas/library/config.h"
#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/gridtools/GridToolsTraits.h"

#if ATLAS_HAVE_ACC
#include "atlas_acc_support/atlas_acc.h"
#endif

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace gridtools {

template <typename gt_DataStore>
struct GridToolsDataStore : ArrayDataStore {
    using Value = typename gt_DataStore::data_t;

    explicit GridToolsDataStore(gt_DataStore const* ds): data_store_(ds) {}

    ~GridToolsDataStore() {
        assert(data_store_);
        delete data_store_;
    }

    void updateDevice() const override {
        assert(data_store_);
        allocateDevice();
        data_store_->clone_to_device();
    }

    void updateHost() const override { data_store_->clone_from_device(); }

    bool valid() const override { return data_store_->valid(); }

    void syncHostDevice() const override { data_store_->sync(); }

    void allocateDevice() const override {
        accMap();
    }

    void deallocateDevice() const override {
        accUnmap();
    }

    bool deviceAllocated() const override { return ATLAS_HAVE_GPU; }

    bool hostNeedsUpdate() const override { return data_store_->host_needs_update(); }

    bool deviceNeedsUpdate() const override { return data_store_->device_needs_update(); }

    void setHostNeedsUpdate(bool v) const override {
        auto state_machine = data_store_->get_storage_ptr()->get_state_machine_ptr_impl();
        if (state_machine) {
            state_machine->m_hnu = v;
        }
    }

    void setDeviceNeedsUpdate(bool v) const override {
        auto state_machine = data_store_->get_storage_ptr()->get_state_machine_ptr_impl();
        if (state_machine) {
            state_machine->m_dnu = v;
        }
    }

    void reactivateDeviceWriteViews() const override { data_store_->reactivate_target_write_views(); }

    void reactivateHostWriteViews() const override { data_store_->reactivate_host_write_views(); }

    void* voidDataStore() override { return static_cast<void*>(const_cast<gt_DataStore*>(data_store_)); }

    void* voidHostData() override { return host_data(); }

    void* voidDeviceData() override { return device_data(); }

    void accMap() const override {
#if ATLAS_HAVE_ACC
        if (not acc_mapped_) {
            ATLAS_ASSERT(deviceAllocated(),"Could not accMap as device data is not allocated");
            atlas_acc_map_data(host_data(), device_data(), bytes());
            acc_mapped_ = true;
        }
#endif
    }

    bool accMapped() const override {
        return acc_mapped_;
    }

    void accUnmap() const override {
#if ATLAS_HAVE_ACC
        if (acc_mapped_) {
            atlas_acc_unmap_data(host_data());
            acc_mapped_ = false;
        }
#endif
    }

private:
    void* host_data() const {
        return ::gridtools::make_host_view<::gridtools::access_mode::read_only>(*data_store_).data();
    }
    void* device_data() const {
#if ATLAS_HAVE_GPU
        return ::gridtools::make_device_view<::gridtools::access_mode::read_only>(*data_store_).data();
#else
        return ::gridtools::make_host_view<::gridtools::access_mode::read_only>(*data_store_).data();
#endif
    }

    size_t bytes() const {
        auto storage_info_ptr = data_store_->get_storage_info_ptr().get();
        return storage_info_ptr->padded_total_length() * sizeof(Value);
    }

private:
    gt_DataStore const* data_store_;
    mutable bool acc_mapped_{false};
};

}  // namespace gridtools
}  // namespace array
}  // namespace atlas
