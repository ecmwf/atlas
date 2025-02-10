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

#include <memory>
#include <vector>

#include "atlas/library/config.h"

#include "atlas/util/Object.h"

#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/DataType.h"
#include "atlas/array_fwd.h"

namespace atlas {
namespace array {

// --------------------------------------------------------------------------------------------
// Forward declarations
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <typename Value>
class ArrayT;
template <typename Value>
class ArrayT_impl;
#endif

// --------------------------------------------------------------------------------------------

class Array : public util::Object {
public:
    Array() = default;
    virtual ~Array();

    static Array* create(array::DataType, const ArrayShape&);

    static Array* create(array::DataType, const ArrayShape&, const ArrayLayout&);

    static Array* create(array::DataType, ArraySpec&&);

    static Array* create(ArraySpec&&);

    virtual size_t footprint() const = 0;

    template <typename Value>
    static Array* create(idx_t size0);
    template <typename Value>
    static Array* create(idx_t size0, idx_t size1);
    template <typename Value>
    static Array* create(idx_t size0, idx_t size1, idx_t size2);
    template <typename Value>
    static Array* create(idx_t size0, idx_t size1, idx_t size2, idx_t size3);
    template <typename Value>
    static Array* create(idx_t size0, idx_t size1, idx_t size2, idx_t size3, idx_t size4);

    template <typename Value>
    static Array* create(const ArrayShape& shape);

    template <typename Value>
    static Array* create(const ArrayShape& shape, const ArrayLayout& layout);

    template <typename Value>
    static Array* wrap(Value* data, const ArrayShape& shape);

    template <typename Value>
    static Array* wrap(Value* data, const ArraySpec& spec);

    idx_t bytes() const { return datatype().size() * spec().allocatedSize(); }

    size_t size() const { return spec_.size(); }

    idx_t rank() const { return spec_.rank(); }

    idx_t stride(idx_t i) const { return spec_.strides()[i]; }

    idx_t shape(idx_t i) const { return spec_.shape()[i]; }

    const ArrayStrides& strides() const { return spec_.strides(); }

    const ArrayStrides& device_strides() const { return spec_.device_strides(); }

    const ArrayShape& shape() const { return spec_.shape(); }

    const std::vector<int>& shapef() const { return spec_.shapef(); }

    const std::vector<int>& stridesf() const { return spec_.stridesf(); }

    const std::vector<int>& device_stridesf() const { return spec_.device_stridesf(); }

    bool contiguous() const { return spec_.contiguous(); }

    bool hasDefaultLayout() const { return spec_.hasDefaultLayout(); }

    virtual array::DataType datatype() const = 0;

    virtual void resize(const ArrayShape& shape) = 0;

    virtual void resize(idx_t size0)                                                     = 0;
    virtual void resize(idx_t size0, idx_t size1)                                        = 0;
    virtual void resize(idx_t size0, idx_t size1, idx_t size2)                           = 0;
    virtual void resize(idx_t size0, idx_t size1, idx_t size2, idx_t size3)              = 0;
    virtual void resize(idx_t size0, idx_t size1, idx_t size2, idx_t size3, idx_t size4) = 0;

    virtual void insert(idx_t idx1, idx_t size1) = 0;

    virtual void dump(std::ostream& os) const = 0;

    virtual void accMap() const = 0;
    virtual void accUnmap() const = 0;
    virtual bool accMapped() const = 0;

    virtual void* storage() { return data_store_->voidDataStore(); }

    virtual const void* storage() const { return data_store_->voidDataStore(); }

    bool valid() const { return data_store_->valid(); }

    void updateDevice() const { data_store_->updateDevice(); }

    void updateHost() const { data_store_->updateHost(); }

    void syncHostDevice() const { data_store_->syncHostDevice(); }

    bool hostNeedsUpdate() const { return data_store_->hostNeedsUpdate(); }

    bool deviceNeedsUpdate() const { return data_store_->deviceNeedsUpdate(); }

    void setHostNeedsUpdate(bool v) const { return data_store_->setHostNeedsUpdate(v); }

    void setDeviceNeedsUpdate(bool v) const { return data_store_->setDeviceNeedsUpdate(v); }

    bool deviceAllocated() const { return data_store_->deviceAllocated(); }

    void allocateDevice() { data_store_->allocateDevice(); }

    void deallocateDevice() { data_store_->deallocateDevice(); }

    void reactivateDeviceWriteViews() const { data_store_->reactivateDeviceWriteViews(); }

    void reactivateHostWriteViews() const { data_store_->reactivateHostWriteViews(); }

    const ArraySpec& spec() const { return spec_; }

    struct CopyPolicy {
        enum class Execution {
            SERIAL=0,
            OMP=1
        };
        bool on_device = false;
        Execution execution {Execution::SERIAL};
    };
    virtual void copy(const Array&, const CopyPolicy&) = 0;
    void copy(const Array& other) { return copy(other,CopyPolicy{}); }

    // -- dangerous methods... You're on your own interpreting the raw data
    template <typename DATATYPE>
    DATATYPE const* host_data() const {
        return data_store_->hostData<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE* host_data() {
        return data_store_->hostData<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE const* device_data() const {
        return data_store_->deviceData<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE* device_data() {
        return data_store_->deviceData<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE const* data() const {
        return data_store_->hostData<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE* data() {
        return data_store_->hostData<DATATYPE>();
    }
    void const* data() const { return data<void>(); }
    void* data() { return data<void>(); }

    const ArrayDataStore& data_store() const { return *data_store_; }

protected:
    Array(ArraySpec&& spec): spec_(std::move(spec)) {}
    ArraySpec spec_;
    std::unique_ptr<ArrayDataStore> data_store_;

    void replace(Array& array) {
        data_store_.swap(array.data_store_);
        spec_ = array.spec_;
    }
};

// --------------------------------------------------------------------------------------------

template <typename Value>
class ArrayT : public Array {
public:
    ArrayT(idx_t size0);
    ArrayT(idx_t size0, idx_t size1);
    ArrayT(idx_t size0, idx_t size1, idx_t size2);
    ArrayT(idx_t size0, idx_t size1, idx_t size2, idx_t size3);
    ArrayT(idx_t size0, idx_t size1, idx_t size2, idx_t size3, idx_t size4);

    ArrayT(ArraySpec&&);

    ArrayT(const ArrayShape&);

    ArrayT(const ArrayShape&, const ArrayAlignment&);

    ArrayT(const ArrayShape&, const ArrayLayout&);

    virtual void insert(idx_t idx1, idx_t size1);

    virtual void resize(const ArrayShape&);

    virtual void resize(idx_t size0);
    virtual void resize(idx_t size0, idx_t size1);
    virtual void resize(idx_t size0, idx_t size1, idx_t size2);
    virtual void resize(idx_t size0, idx_t size1, idx_t size2, idx_t size3);
    virtual void resize(idx_t size0, idx_t size1, idx_t size2, idx_t size3, idx_t size4);

    virtual void copy(const Array&, const CopyPolicy&);

    virtual array::DataType datatype() const { return array::DataType::create<Value>(); }

    virtual void dump(std::ostream& os) const;

    // This constructor is used through the Array::create() or the Array::wrap()
    // methods
    ArrayT(ArrayDataStore*, const ArraySpec&);

    virtual size_t footprint() const;

    virtual void accMap() const;
    virtual void accUnmap() const;
    virtual bool accMapped() const;

    using Array::host_data;
    using Array::device_data;
    using Array::data;

    Value const* host_data() const {
        return data_store_->hostData<Value>();
    }
    Value* host_data() {
        return data_store_->hostData<Value>();
    }
    Value const* device_data() const {
        return data_store_->deviceData<Value>();
    }
    Value* device_data() {
        return data_store_->deviceData<Value>();
    }
    Value const* data() const {
        return data_store_->hostData<Value>();
    }
    Value* data() {
        return data_store_->hostData<Value>();
    }
private:
    template <typename T>
    friend class ArrayT_impl;
};

extern template class ArrayT<float>;
extern template class ArrayT<double>;
extern template class ArrayT<int>;
extern template class ArrayT<long>;

}  // namespace array
}  // namespace atlas
