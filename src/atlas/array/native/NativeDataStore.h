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

#include <algorithm>  // std::fill
#include <limits>   // std::numeric_limits<T>::signaling_NaN
#include <string>
#include <iostream>

#include "pluto/pluto.h"

#include "eckit/utils/Translator.h"

#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/ArraySpec.h"
#include "atlas/library.h"
#include "atlas/library/config.h"
#include "atlas/parallel/acc/acc.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"

#include "atlas/runtime/Log.h"


#define ATLAS_ACC_DEBUG 0

namespace {
bool getEnv(const std::string& env, bool default_value) {
    const char* cenv = ::getenv(env.c_str());
    if (cenv != nullptr) {
        return eckit::Translator<std::string, bool>()(cenv);
    }
    return default_value;
}
#if ATLAS_HAVE_GPU
bool DEVICE_MEMORY() {
    static bool v = getEnv("ATLAS_SIMULATE_DEVICE_MEMORY", true) || pluto::devices() > 0;
    return v;
}
#else
bool DEVICE_MEMORY() {
  static bool ATLAS_ALLOCATE_DEVICE = getEnv("ATLAS_SIMULATE_DEVICE_MEMORY",false);
  return ATLAS_ALLOCATE_DEVICE;
}
}
#endif

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace native {

template <typename Value>
static constexpr Value invalid_value() {
    return std::numeric_limits<Value>::has_signaling_NaN ? std::numeric_limits<Value>::signaling_NaN()
           : std::numeric_limits<Value>::has_quiet_NaN   ? std::numeric_limits<Value>::quiet_NaN()
           : std::numeric_limits<Value>::has_infinity    ? std::numeric_limits<Value>::infinity()
                                                         : std::numeric_limits<Value>::max();
}

#if ATLAS_INIT_SNAN
template <typename Value>
void initialise(Value array[], size_t size) {
    std::fill_n(array, size, invalid_value<Value>());
}
#else
template <typename Value>
void initialise(Value[], size_t) {}
#endif

static bool test_device_is_shared() {
    // This is currently a bit hacky, to be improved!
    auto& device_resource = *pluto::device::get_default_resource();

    if (device_resource == *pluto::pinned_resource()) {
        return true;
    }
    if (device_resource == *pluto::managed_resource()) {
        return true;
    }
    if (device_resource == *pluto::pinned_pool_resource()) {
        return true;
    }
    if (device_resource == *pluto::managed_pool_resource()) {
        return true;
    }

    // If above shortcuts failed, we can do an introspection on a dummy allocation
    if (DEVICE_MEMORY()) {
        void* device_ptr = device_resource.allocate(8);
        bool is_shared = pluto::is_managed(device_ptr) || (pluto::is_pinned(device_ptr));
        device_resource.deallocate(device_ptr,8);
        return is_shared;
    }
    return false;
}

static pluto::memory_resource* host_resource(bool shared) {
    if (shared) {
        return pluto::device::get_default_resource();
    }
    else {
        return pluto::host::get_default_resource();
    }
}

static pluto::memory_resource* device_resource(bool shared) {
    if (shared) {
        return pluto::null_memory_resource();
    }
    else {
        return pluto::device::get_default_resource();
    }
}

static pluto::memory_resource* device_mapped_host_resource() {
    if (pluto::has_registered_resource("atlas::array::device_mapped_host_resource")) {
        return pluto::get_registered_resource("atlas::array::device_mapped_host_resource");
    }
    return pluto::host::get_default_resource();
}

template <typename Value>
class DataStore final : public ArrayDataStore {
public:
    DataStore(size_t size): size_(size),
        is_shared_data_{test_device_is_shared()},
        host_allocator_{host_resource(is_shared_data_)},
        device_allocator_{device_resource(is_shared_data_)} {

        label_ = array::label::get();
        allocateHost();
        initialise(host_data_, size_);
        if (DEVICE_MEMORY()) {
            device_updated_ = false;
        }
        else {
            device_data_ = host_data_;
            device_mapped_host_data_ = host_data_;
        }
    }

    virtual ~DataStore() {
        deallocateDevice();
        deallocateHost();
    }

    void updateDevice() const override {
        if (DEVICE_MEMORY()) {
            if (not device_allocated_) {
                allocateDevice();
            }
#if ATLAS_HAVE_TRACE
            std::unique_ptr<atlas::Trace> trace;
            if (atlas::Library::instance().traceMemory()) {
                trace.reset(new atlas::Trace(Here(), "update_device", {"memory","memory.transfer", "memory.transfer.h2d"}));
            }
#endif
            pluto::copy_host_to_device(label_, device_data_, host_data_, size_);
        }
        device_updated_ = true;
    }

    void updateHost() const override {
        ATLAS_DEBUG_VAR(DEVICE_MEMORY());
        if (DEVICE_MEMORY()) {
            ATLAS_DEBUG_VAR(device_allocated_);
            if (device_allocated_) {
#if ATLAS_HAVE_TRACE
                std::unique_ptr<atlas::Trace> trace;
                if (atlas::Library::instance().traceMemory()) {
                    trace.reset(new atlas::Trace(Here(), "update_host", {"memory","memory.transfer", "memory.transfer.d2h"}));
                }
#endif
                ATLAS_DEBUG("copy_device_to_host");
                pluto::copy_device_to_host(label_, host_data_, device_data_, size_);
            }
        }
        host_updated_ = true;
    }

    bool valid() const override { return true; }

    void syncHostDevice() const override {
        if (host_updated_ and device_updated_) {
            return; // nothing to do
        }
        if (not (host_updated_ or device_updated_)) {
            throw_AssertionFailed("syncHostDevice() could not figure out which of host or device is up to date. "
                                  "Probably it was forgotten to use setDeviceNeedsUpdate(true) or setDeviceNeedsUpdate(true)",
                                  Here());
        }

        if (not device_updated_) {
            updateDevice();
        }
        else if (not host_updated_) {
            updateHost();
        }
    }

    void syncHost() const override {
        if (not host_updated_) {
            ATLAS_ASSERT(device_updated_,"Unexpected state when calling syncHost() because (deviceNeedsUpdate() and hostNeedsUpdate())");
            updateHost();
        }
    }

    void syncDevice() const override {
        if (not device_updated_) {
            ATLAS_ASSERT(host_updated_,"Unexpected state when calling syncDevice() because (deviceNeedsUpdate() and hostNeedsUpdate())");
            if (not host_updated_) {
                throw_AssertionFailed("Unexpected state when calling syncDevice() because (deviceNeedsUpdate() and hostNeedsUpdate())",
                Here());
            }
            updateDevice();
        }
    }

    bool deviceAllocated() const override { return device_allocated_; }

    void allocateDevice() const override {
        if (DEVICE_MEMORY()) {
            if (device_allocated_) {
                return;
            }
            if (size_) {
                if(is_shared_data_) {
                    device_data_ = pluto::get_registered_device_pointer(host_data_);
                }
                else {
                    #if ATLAS_HAVE_TRACE
                        std::unique_ptr<atlas::Trace> trace;
                        if (atlas::Library::instance().traceMemory()) {
                            trace.reset(new atlas::Trace(Here(), "allocate_device", {"memory","memory.allocator","memory.allocator.device","memory.allocator.device.allocate"}));
                        }
                    #endif
                    device_data_ = device_allocator_.allocate(label_, size_);
                }
                if (pluto::devices() > 0) {
                    device_mapped_host_data_ = host_data_;
                    ATLAS_ASSERT(pluto::is_device_accessible(device_data_));
                }
                else {
                    device_mapped_host_data_ = device_data_;
                }
                device_allocated_ = true;
                accMap();
            }
        }
    }

    void deallocateDevice() const override {
        if (device_allocated_) {
            accUnmap();
            if (not is_shared_data_) {
#if ATLAS_HAVE_TRACE
                std::unique_ptr<atlas::Trace> trace;
                if (atlas::Library::instance().traceMemory()) {
                    trace.reset(new atlas::Trace(Here(), "deallocate_device", {"memory","memory.allocator","memory.allocator.device","memory.allocator.device.deallocate"}));
                }
#endif
                device_allocator_.deallocate(label_, device_data_,size_);
            }
            device_mapped_host_data_ = nullptr;
            device_data_ = nullptr;
            device_allocated_ = false;
            device_updated_ = false;
        }
    }

    bool hostNeedsUpdate() const override { return (not host_updated_); }

    bool deviceNeedsUpdate() const override { return (not device_updated_); }

    void setHostNeedsUpdate(bool v) const override { host_updated_ = (not v); }

    void setDeviceNeedsUpdate(bool v) const override { device_updated_ = (not v); }

    void reactivateDeviceWriteViews() const override {}

    void reactivateHostWriteViews() const override {}

    void* voidDataStore() override { return static_cast<void*>(host_data_); }

    void* voidHostData() override { return static_cast<void*>(host_data_); }

    void* voidDeviceData() override { return static_cast<void*>(device_data_); }

    void* voidDeviceMappedHostData() override { return static_cast<void*>(device_mapped_host_data_); }

    void accMap() const override {
        if ((not acc_mapped_ && acc::devices())) {
            if (is_shared_data_) {
                // nvidia compiler should not accMap for managed memory
                if (pluto::is_managed(host_data_) && acc::compiler_id() == acc::CompilerId::nvidia) {
                  return;
                }
            }
            ATLAS_ASSERT(deviceAllocated(),"Could not accMap as device data is not allocated");
            ATLAS_ASSERT(!atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            if constexpr(ATLAS_ACC_DEBUG) {
                std::cout << "               + acc_map_data(hostptr:"<<host_data_<<", device:"<<device_data_<<", bytes:"<<footprint()<<")" <<std::endl;
            }
            atlas::acc::map((void*)host_data_, (void*)device_data_, size_ * sizeof(Value));
            ATLAS_ASSERT(atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            ATLAS_ASSERT(atlas::acc::deviceptr(host_data_) == device_data_);
            acc_mapped_ = true;
        }
    }

    bool accMapped() const override {
        return acc_mapped_;
    }

    void accUnmap() const override {
        if (acc_mapped_) {
            if constexpr(ATLAS_ACC_DEBUG) {
                std::cout << "               - acc_unmap_data(hostptr:"<<host_data_<<", device:"<<device_data_<<", bytes:"<<footprint()<<")" <<std::endl;
            }
            ATLAS_ASSERT(atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            atlas::acc::unmap(host_data_);
            ATLAS_ASSERT(!atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            acc_mapped_ = false;
        }
    }


private:

    void allocateHost() {
        if (size_ > 0) {
            #if ATLAS_HAVE_TRACE
                std::unique_ptr<atlas::Trace> trace;
                if (atlas::Library::instance().traceMemory()) {
                    trace.reset(new atlas::Trace(Here(), "allocate_host", {"memory","memory.allocator","memory.allocator.host","memory.allocator.host.allocate"}));
                }
            #endif
            host_data_ = host_allocator_.allocate(label_, size_);
        }
        else {
            host_data_ = nullptr;
        }
    }

    void deallocateHost() {
        if (host_data_) {
            #if ATLAS_HAVE_TRACE
                std::unique_ptr<atlas::Trace> trace;
                if (atlas::Library::instance().traceMemory()) {
                    trace.reset(new atlas::Trace(Here(), "deallocate_host", {"memory","memory.allocator","memory.allocator.host","memory.allocator.host.deallocate"}));
                }
            #endif
            host_allocator_.deallocate(label_, host_data_, size_);
            host_data_ = nullptr;
        }
    }

    size_t footprint() const { return sizeof(Value) * size_; }

    size_t size_;
    Value* host_data_;
    mutable Value* device_data_{nullptr};
    mutable Value* device_mapped_host_data_{nullptr};

    mutable bool host_updated_{true};
    mutable bool device_updated_{true};
    mutable bool device_allocated_{false};
    mutable bool acc_mapped_{false};
    bool is_shared_data_{false};

    pluto::allocator<Value> host_allocator_;
    mutable pluto::allocator<Value> device_allocator_;

    std::string label_;
};

//------------------------------------------------------------------------------

template <typename Value>
class WrappedDataStore final : public ArrayDataStore {
public:

    void init_device(bool contiguous) {
        if (DEVICE_MEMORY()) {
            device_updated_ = false;
        }
        else {
            device_data_ = host_data_;
            device_mapped_host_data_ = host_data_;
        }
    }

    WrappedDataStore(Value* host_data, size_t size): host_data_(host_data), size_(size),
        device_allocator_{pluto::device::get_default_resource()},
        device_mapped_host_allocator_{device_mapped_host_resource()} {
        label_ = array::label::get();
        init_device(contiguous_);
    }

    WrappedDataStore(Value* host_data, const ArraySpec& spec):
        host_data_(host_data),
        size_(spec.size()),
        device_allocator_{pluto::device::get_default_resource()},
        device_mapped_host_allocator_{device_mapped_host_resource()} {
        label_ = array::label::get();
        contiguous_ = spec.contiguous();
        init_device(contiguous_);
        if (! contiguous_) {
            int break_idx = 0;
            size_t shp_mult_rhs = spec.shape()[spec.rank() - 1];
            for (int i = spec.rank() - 1; i > 0; i--) {
                if (shp_mult_rhs != spec.strides()[i - 1]) {
                    break_idx = i - 1;
                    break;
                }
                shp_mult_rhs *= spec.shape()[i - 1];
            }
            size_t shp_mult_lhs = spec.shape()[0];
            for (int i = 1; i <= break_idx; i++) {
                shp_mult_lhs *= spec.shape()[i];
            }
            if (spec.strides()[spec.rank() - 1] > 1) {
                memcpy_h2d_pitch_ = 1;
                memcpy_d2h_pitch_ = spec.strides()[spec.rank() - 1];
                memcpy_width_ = 1;
                memcpy_height_ = spec.shape()[0] * spec.device_strides()[0];
            }
            else {
                memcpy_h2d_pitch_ = shp_mult_rhs;
                memcpy_d2h_pitch_ = spec.strides()[break_idx];
                memcpy_width_ = shp_mult_rhs;
                memcpy_height_ = shp_mult_lhs;
            }
        }
    }

    virtual ~WrappedDataStore() {
        deallocateDevice();
    }

    void updateDevice() const override {
        allocateDevice();
        if (DEVICE_MEMORY()) {
#if ATLAS_HAVE_TRACE
            std::unique_ptr<atlas::Trace> trace;
            if (atlas::Library::instance().traceMemory()) {
                trace.reset(new atlas::Trace(Here(), "update_device [host=wrapped]", {"memory","memory.transfer", "memory.transfer.h2d"}));
            }
#endif
            if (contiguous_) {
                pluto::copy_host_to_device(label_, device_data_, host_data_, size_);
            }
            else {
                pluto::copy_host_to_device_2D(
                    label_,
                    device_data_, memcpy_h2d_pitch_,
                    host_data_, memcpy_d2h_pitch_,
                    memcpy_width_, memcpy_height_);
            }
        }
        device_updated_ = true;
    }

    void updateHost() const override {
        if (DEVICE_MEMORY()) {
            if (device_allocated_) {
#if ATLAS_HAVE_TRACE
                std::unique_ptr<atlas::Trace> trace;
                if (atlas::Library::instance().traceMemory()) {
                    trace.reset(new atlas::Trace(Here(), "update_host [host=wrapped]", {"memory","memory.transfer", "memory.transfer.d2h"}));
                }
#endif
                if (contiguous_) {
                    pluto::copy_device_to_host(label_, host_data_, device_data_, size_);
                }
                else {
                    pluto::copy_device_to_host_2D(
                        label_,
                        host_data_, memcpy_d2h_pitch_ ,
                        device_data_, memcpy_h2d_pitch_,
                        memcpy_width_, memcpy_height_);
                }
            }
        }
        host_updated_ = true;
    }


    bool valid() const override { return true; }

    void syncHostDevice() const override {
        if (host_updated_ and device_updated_) {
            return; // nothing to do
        }
        if (not (host_updated_ or device_updated_)) {
            throw_AssertionFailed("syncHostDevice() could not figure out which of host or device is up to date. "
                                  "Probably it was forgotten to use setDeviceNeedsUpdate(true) or setDeviceNeedsUpdate(true)",
                                  Here());
        }

        if (not device_updated_) {
            updateDevice();
        }
        else if (not host_updated_) {
            updateHost();
        }
    }

    void syncHost() const override {
        if (not host_updated_) {
            ATLAS_ASSERT(device_updated_,"Unexpected state when calling syncHost() because (deviceNeedsUpdate() and hostNeedsUpdate())");
            updateHost();
        }
    }

    void syncDevice() const override {
        if (not device_updated_) {
            ATLAS_ASSERT(host_updated_,"Unexpected state when calling syncDevice() because (deviceNeedsUpdate() and hostNeedsUpdate())");
            updateDevice();
        }
    }

    bool deviceAllocated() const override { return device_allocated_; }

    void allocateDevice() const override {
        if (DEVICE_MEMORY()) {
            if (device_allocated_) {
                return;
            }
            if (size_) {
#if ATLAS_HAVE_TRACE
                std::unique_ptr<atlas::Trace> trace;
                if (atlas::Library::instance().traceMemory()) {
                    trace.reset(new atlas::Trace(Here(), "allocate_device [host=wrapped]", {"memory","memory.allocator","memory.allocator.device","memory.allocator.device.allocate"}));
                }
#endif
                device_data_ = device_allocator_.allocate(label_, size_);
                device_allocated_ = true;
                if (pluto::devices() == 0) {
                    device_mapped_host_data_ = device_data_;
                }
                else if (contiguous_) {
                    device_mapped_host_data_ = host_data_;
                    accMap();
                }
                else {
                    device_mapped_host_data_ = nullptr;
                }
            }
        }
    }

    void allocateDeviceMappedHostData() const {
        if (not device_mapped_host_data_) {
            device_mapped_host_data_ = device_mapped_host_allocator_.allocate("device_mapped_host("+label_+")", size_);
        }
    }

    void deallocateDeviceMappedHostData() const {
        if (device_mapped_host_data_) {
            if (device_mapped_host_data_ != device_data_ && device_mapped_host_data_ != host_data_) {
                device_mapped_host_allocator_.deallocate("device_mapped_host("+label_+")",device_mapped_host_data_, size_);
                device_mapped_host_data_ = nullptr;
            }
        }
    }

    void deallocateDevice() const override {
        if (device_allocated_) {
            accUnmap();
#if ATLAS_HAVE_TRACE
            std::unique_ptr<atlas::Trace> trace;
            if (atlas::Library::instance().traceMemory()) {
                trace.reset(new atlas::Trace(Here(), "deallocate_device [host=wrapped]", {"memory","memory.allocator","memory.allocator.device","memory.allocator.device.deallocate"}));
            }
#endif
            deallocateDeviceMappedHostData();

            device_allocator_.deallocate(label_, device_data_, size_);
            device_data_ = nullptr;
            device_allocated_ = false;
            device_updated_ = false;
        }
    }

    bool hostNeedsUpdate() const override { return (not host_updated_); }

    bool deviceNeedsUpdate() const override { return (not device_updated_); }

    void setHostNeedsUpdate(bool v) const override { host_updated_ = (not v); }

    void setDeviceNeedsUpdate(bool v) const override { device_updated_ = (not v); }

    void reactivateDeviceWriteViews() const override {}

    void reactivateHostWriteViews() const override {}

    void* voidDataStore() override { return static_cast<void*>(host_data_); }

    void* voidHostData() override { return static_cast<void*>(host_data_); }

    void* voidDeviceData() override { return static_cast<void*>(device_data_); }

    void* voidDeviceMappedHostData() override {
        if (device_mapped_host_data_) {
            return static_cast<void*>(device_mapped_host_data_);
        }
        allocateDeviceMappedHostData();
        accMap();
        return static_cast<void*>(device_mapped_host_data_);
    }

    void accMap() const override {
        if (not acc_mapped_ && acc::devices()) {
            auto do_acc_map = [&](void* hostptr, void* deviceptr) {
                ATLAS_ASSERT(deviceAllocated(),"Could not accMap as device data is not allocated");
                ATLAS_ASSERT(!atlas::acc::is_present(hostptr, size_ * sizeof(Value)));
                if constexpr(ATLAS_ACC_DEBUG) {
                    std::cout << "               + acc_map_data(hostptr:"<<hostptr<<", deviceptr:"<<deviceptr<<", bytes:"<<size_ * sizeof(Value)<<")" <<std::endl;
                }
                atlas::acc::map((void*)hostptr, (void*)deviceptr, size_ * sizeof(Value));
                ATLAS_ASSERT(atlas::acc::is_present(hostptr, size_ * sizeof(Value)));
                ATLAS_ASSERT(atlas::acc::deviceptr(hostptr) == deviceptr);
                acc_mapped_ = true;
            };
            if (contiguous_) {
                do_acc_map(host_data_, device_data_);
            }
            else {
                do_acc_map(device_mapped_host_data_, device_data_);
            }
        }
    }

    bool accMapped() const override {
        return acc_mapped_;
    }

    void accUnmap() const override {
        if (acc_mapped_) {
            auto do_acc_unmap = [&](void* hostptr, void* deviceptr) {
                ATLAS_ASSERT(atlas::acc::is_present(hostptr, size_ * sizeof(Value)));
                if constexpr(ATLAS_ACC_DEBUG) {
                    std::cout << "               - acc_unmap_data(hostptr:"<<hostptr<<", deviceptr:"<<deviceptr<<", bytes:"<<size_ * sizeof(Value)<<")" <<std::endl;
                }
                atlas::acc::unmap(hostptr);
                if constexpr(ATLAS_ACC_DEBUG) {
                    std::cout << "               - acc_unmap_data(hostptr:"<<hostptr<<", deviceptr:"<<deviceptr<<", bytes:"<<size_ * sizeof(Value)<<") ... done" <<std::endl;
                }
                acc_mapped_ = false;
            };
            if (contiguous_) {
                do_acc_unmap(host_data_, device_data_);
            }
            else if(device_mapped_host_data_) {
                do_acc_unmap(device_mapped_host_data_, device_data_);
            }
        }
    }

private:
    Value* host_data_;
    size_t size_;
    mutable Value* device_data_;
    mutable Value* device_mapped_host_data_{nullptr};

    bool contiguous_{true};
    size_t memcpy_h2d_pitch_;
    size_t memcpy_d2h_pitch_;
    size_t memcpy_height_;
    size_t memcpy_width_;

    mutable bool host_updated_{true};
    mutable bool device_updated_{true};
    mutable bool device_allocated_{false};
    mutable bool acc_mapped_{false};

    mutable pluto::allocator<Value> device_allocator_;
    mutable pluto::allocator<Value> device_mapped_host_allocator_;

    std::string label_;
};

}  // namespace native
}  // namespace array
}  // namespace atlas
