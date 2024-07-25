/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>

#include "eckit/config/Parametrisation.h"

#include "atlas/array.h"
#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/helpers/ArrayInitializer.h"
#include "atlas/array/helpers/ArrayWriter.h"
#include "atlas/array/native/NativeDataStore.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Exception.h"

using namespace atlas::array::helpers;

namespace atlas {
namespace array {

template <typename Value>
Array* Array::create(idx_t dim0, const eckit::Parametrisation& param) {
    return new ArrayT<Value>(dim0, param);
}
template <typename Value>
Array* Array::create(idx_t dim0, idx_t dim1, const eckit::Parametrisation& param) {
    return new ArrayT<Value>(dim0, dim1, param);
}
template <typename Value>
Array* Array::create(idx_t dim0, idx_t dim1, idx_t dim2, const eckit::Parametrisation& param) {
    return new ArrayT<Value>(dim0, dim1, dim2, param);
}
template <typename Value>
Array* Array::create(idx_t dim0, idx_t dim1, idx_t dim2, idx_t dim3, const eckit::Parametrisation& param) {
    return new ArrayT<Value>(dim0, dim1, dim2, dim3, param);
}
template <typename Value>
Array* Array::create(idx_t dim0, idx_t dim1, idx_t dim2, idx_t dim3, idx_t dim4,
        const eckit::Parametrisation& param) {
    return new ArrayT<Value>(dim0, dim1, dim2, dim3, dim4, param);
}
template <typename Value>
Array* Array::create(const ArrayShape& shape, const eckit::Parametrisation& param) {
    return new ArrayT<Value>(shape, param);
}
template <typename Value>
Array* Array::create(const ArrayShape& shape, const ArrayLayout& layout,
        const eckit::Parametrisation& param) {
    return new ArrayT<Value>(shape, layout, param);
}
template <typename Value>
Array* Array::wrap(Value* data, const ArrayShape& shape, const eckit::Parametrisation& param) {
    size_t size = 1;
    for (int i = 0; i < shape.size(); ++i) {
        size *= shape[i];
    }
    return new ArrayT<Value>(new native::WrappedDataStore<Value>(data, size, param), shape, param);
}
template <typename Value>
Array* Array::wrap(Value* data, const ArraySpec& spec, const eckit::Parametrisation& param) {
    size_t size = spec.size();
    return new ArrayT<Value>(new native::WrappedDataStore<Value>(data, spec, param), spec, param);
}

Array::~Array() = default;

Array* Array::create(DataType datatype, const ArrayShape& shape, const eckit::Parametrisation& param) {
    switch (datatype.kind()) {
        case DataType::KIND_REAL64:
            return new ArrayT<double>(shape, param);
        case DataType::KIND_REAL32:
            return new ArrayT<float>(shape, param);
        case DataType::KIND_INT32:
            return new ArrayT<int>(shape, param);
        case DataType::KIND_INT64:
            return new ArrayT<long>(shape, param);
        case DataType::KIND_UINT64:
            return new ArrayT<unsigned long>(shape, param);
        default: {
            std::stringstream err;
            err << "data kind " << datatype.kind() << " not recognised.";
            throw_NotImplemented(err.str(), Here());
        }
    }
}

Array* Array::create(DataType datatype, ArraySpec&& spec, const eckit::Parametrisation& param) {
    switch (datatype.kind()) {
        case DataType::KIND_REAL64:
            return new ArrayT<double>(std::move(spec), param);
        case DataType::KIND_REAL32:
            return new ArrayT<float>(std::move(spec), param);
        case DataType::KIND_INT32:
            return new ArrayT<int>(std::move(spec), param);
        case DataType::KIND_INT64:
            return new ArrayT<long>(std::move(spec), param);
        case DataType::KIND_UINT64:
            return new ArrayT<unsigned long>(std::move(spec), param);
        default: {
            std::stringstream err;
            err << "data kind " << datatype.kind() << " not recognised.";
            throw_NotImplemented(err.str(), Here());
        }
    }
}

Array* Array::create(ArraySpec&& spec, const eckit::Parametrisation& param) {
    return create(spec.datatype(), std::move(spec), param);
}


template <typename Value>
ArrayT<Value>::ArrayT(ArrayDataStore* ds, const ArraySpec& spec, const eckit::Parametrisation& param) {
    data_store_ = std::unique_ptr<ArrayDataStore>(ds);
    spec_       = spec;
    param.get("host_memory_mapped",mapped_);
}

template <typename Value>
ArrayT<Value>::ArrayT(idx_t dim0, const eckit::Parametrisation& param) {
    spec_       = ArraySpec(make_shape(dim0));
    data_store_ = std::make_unique<native::DataStore<Value>>(spec_.size(), param);
    param.get("host_memory_mapped",mapped_);
}
template <typename Value>
ArrayT<Value>::ArrayT(idx_t dim0, idx_t dim1, const eckit::Parametrisation& param) {
    spec_       = ArraySpec(make_shape(dim0, dim1));
    data_store_ = std::make_unique<native::DataStore<Value>>(spec_.size(), param);
    param.get("host_memory_mapped",mapped_);
}
template <typename Value>
ArrayT<Value>::ArrayT(idx_t dim0, idx_t dim1, idx_t dim2, const eckit::Parametrisation& param) {
    spec_       = ArraySpec(make_shape(dim0, dim1, dim2));
    data_store_ = std::make_unique<native::DataStore<Value>>(spec_.size(), param);
    param.get("host_memory_mapped",mapped_);
}
template <typename Value>
ArrayT<Value>::ArrayT(idx_t dim0, idx_t dim1, idx_t dim2, idx_t dim3, const eckit::Parametrisation& param) {
    spec_       = ArraySpec(make_shape(dim0, dim1, dim2, dim3));
    data_store_ = std::make_unique<native::DataStore<Value>>(spec_.size(), param);
    param.get("host_memory_mapped",mapped_);
}
template <typename Value>
ArrayT<Value>::ArrayT(idx_t dim0, idx_t dim1, idx_t dim2, idx_t dim3, idx_t dim4,
        const eckit::Parametrisation& param) {
    spec_       = ArraySpec(make_shape(dim0, dim1, dim2, dim3, dim4));
    data_store_ = std::make_unique<native::DataStore<Value>>(spec_.size(), param);
    param.get("host_memory_mapped",mapped_);
}

template <typename Value>
ArrayT<Value>::ArrayT(const ArrayShape& shape, const eckit::Parametrisation& param) {
    ATLAS_ASSERT(shape.size() > 0);
    size_t size = 1;
    for (size_t j = 0; j < shape.size(); ++j) {
        size *= size_t(shape[j]);
    }
    data_store_ = std::make_unique<native::DataStore<Value>>(size, param);
    spec_       = ArraySpec(shape);
    param.get("host_memory_mapped",mapped_);
}

template <typename Value>
ArrayT<Value>::ArrayT(const ArrayShape& shape, const ArrayAlignment& alignment,
        const eckit::Parametrisation& param) {
    spec_       = ArraySpec(shape, alignment);
    data_store_ = std::make_unique<native::DataStore<Value>>(spec_.allocatedSize(), param);
    param.get("host_memory_mapped",mapped_);
}

template <typename Value>
ArrayT<Value>::ArrayT(const ArrayShape& shape, const ArrayLayout& layout, 
        const eckit::Parametrisation& param) {
    spec_       = ArraySpec(shape);
    data_store_ = std::make_unique<native::DataStore<Value>>(spec_.size(), param);
    for (size_t j = 0; j < layout.size(); ++j) {
        ATLAS_ASSERT(spec_.layout()[j] == layout[j]);
    }
    param.get("host_memory_mapped",mapped_);
}

template <typename Value>
ArrayT<Value>::ArrayT(ArraySpec&& spec, const eckit::Parametrisation& param): Array(std::move(spec)) {
    data_store_ = std::make_unique<native::DataStore<Value>>(spec_.allocatedSize(), param);
    param.get("host_memory_mapped",mapped_);
}

template <typename Value>
void ArrayT<Value>::resize(const ArrayShape& _shape) {
    if (rank() != static_cast<idx_t>(_shape.size())) {
        std::stringstream msg;
        msg << "Cannot resize existing Array with rank " << rank() << " with a shape of rank " << _shape.size();
        throw_Exception(msg.str(), Here());
    }

    Array* resized = Array::create<Value>(_shape);

    array_initializer::apply(*this,*resized);
    
    replace(*resized);
    delete resized;
}


template <typename Value>
void ArrayT<Value>::copy(const Array& other, const CopyPolicy&) {
    array_initializer::apply(other,*this);
}

template <typename Value>
void ArrayT<Value>::insert(idx_t idx1, idx_t size1) {
    ArrayShape nshape = shape();
    if (idx1 > nshape[0]) {
        throw_Exception("Cannot insert into an array at a position beyond its size", Here());
    }
    nshape[0] += size1;

    Array* resized = Array::create<Value>(nshape);

    array_initializer_partitioned<0>::apply(*this, *resized, idx1, size1);
    replace(*resized);
    delete resized;
}

template <typename Value>
void ArrayT<Value>::resize(idx_t size1) {
    resize(make_shape(size1));
}

template <typename Value>
void ArrayT<Value>::resize(idx_t size1, idx_t size2) {
    resize(make_shape(size1, size2));
}

template <typename Value>
void ArrayT<Value>::resize(idx_t size1, idx_t size2, idx_t size3) {
    resize(make_shape(size1, size2, size3));
}

template <typename Value>
void ArrayT<Value>::resize(idx_t size1, idx_t size2, idx_t size3, idx_t size4) {
    resize(make_shape(size1, size2, size3, size4));
}

template <typename Value>
void ArrayT<Value>::resize(idx_t size1, idx_t size2, idx_t size3, idx_t size4, idx_t size5) {
    resize(make_shape(size1, size2, size3, size4, size5));
}

template <typename Value>
void ArrayT<Value>::dump(std::ostream& out) const {
    switch (rank()) {
        case 1:
            make_host_view<const Value, 1>(*this).dump(out);
            break;
        case 2:
            make_host_view<const Value, 2>(*this).dump(out);
            break;
        case 3:
            make_host_view<const Value, 3>(*this).dump(out);
            break;
        case 4:
            make_host_view<const Value, 4>(*this).dump(out);
            break;
        case 5:
            make_host_view<const Value, 5>(*this).dump(out);
            break;
        case 6:
            make_host_view<const Value, 6>(*this).dump(out);
            break;
        case 7:
            make_host_view<const Value, 7>(*this).dump(out);
            break;
        case 8:
            make_host_view<const Value, 8>(*this).dump(out);
            break;
        case 9:
            make_host_view<const Value, 9>(*this).dump(out);
            break;
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

//------------------------------------------------------------------------------

template <typename Value>
size_t ArrayT<Value>::footprint() const {
    size_t size = sizeof(*this);
    size += bytes();
    if (not contiguous()) {
        ATLAS_NOTIMPLEMENTED;
    }
    return size;
}

template <typename Value>
void ArrayT<Value>::accMap() const {
    data_store_->accMap();
}

template <typename Value>
void ArrayT<Value>::accUnmap() const {
    data_store_->accUnmap();
}

template <typename Value>
bool ArrayT<Value>::accMapped() const {
    return data_store_->accMapped();
}

//------------------------------------------------------------------------------

template Array* Array::create<int>(idx_t, const eckit::Parametrisation&);
template Array* Array::create<long>(idx_t, const eckit::Parametrisation&);
template Array* Array::create<float>(idx_t, const eckit::Parametrisation&);
template Array* Array::create<double>(idx_t, const eckit::Parametrisation&);
template Array* Array::create<long unsigned>(idx_t, const eckit::Parametrisation&);

template Array* Array::create<int>(idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<long>(idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<float>(idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<double>(idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<long unsigned>(idx_t, idx_t, const eckit::Parametrisation&);

template Array* Array::create<int>(idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<long>(idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<float>(idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<double>(idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<long unsigned>(idx_t, idx_t, idx_t, const eckit::Parametrisation&);

template Array* Array::create<int>(idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<long>(idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<float>(idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<double>(idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<long unsigned>(idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);

template Array* Array::create<int>(idx_t, idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<long>(idx_t, idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<float>(idx_t, idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<double>(idx_t, idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);
template Array* Array::create<long unsigned>(idx_t, idx_t, idx_t, idx_t, idx_t, const eckit::Parametrisation&);

template Array* Array::create<int>(const ArrayShape&, const eckit::Parametrisation&);
template Array* Array::create<long>(const ArrayShape&, const eckit::Parametrisation&);
template Array* Array::create<float>(const ArrayShape&, const eckit::Parametrisation&);
template Array* Array::create<double>(const ArrayShape&, const eckit::Parametrisation&);
template Array* Array::create<long unsigned>(const ArrayShape&, const eckit::Parametrisation&);

template Array* Array::create<int>(const ArrayShape&, const ArrayLayout&, const eckit::Parametrisation&);
template Array* Array::create<long>(const ArrayShape&, const ArrayLayout&, const eckit::Parametrisation&);
template Array* Array::create<float>(const ArrayShape&, const ArrayLayout&, const eckit::Parametrisation&);
template Array* Array::create<double>(const ArrayShape&, const ArrayLayout&, const eckit::Parametrisation&);
template Array* Array::create<long unsigned>(const ArrayShape&, const ArrayLayout&, const eckit::Parametrisation&);

template Array* Array::wrap<int>(int*, const ArrayShape&, const eckit::Parametrisation&);
template Array* Array::wrap<long>(long*, const ArrayShape&, const eckit::Parametrisation&);
template Array* Array::wrap<float>(float*, const ArrayShape&, const eckit::Parametrisation&);
template Array* Array::wrap<double>(double*, const ArrayShape&, const eckit::Parametrisation&);
template Array* Array::wrap<long unsigned>(long unsigned*, const ArrayShape&, const eckit::Parametrisation&);

template Array* Array::wrap<int>(int*, const ArraySpec&, const eckit::Parametrisation&);
template Array* Array::wrap<long>(long*, const ArraySpec&, const eckit::Parametrisation&);
template Array* Array::wrap<float>(float*, const ArraySpec&, const eckit::Parametrisation&);
template Array* Array::wrap<double>(double*, const ArraySpec&, const eckit::Parametrisation&);
template Array* Array::wrap<long unsigned>(long unsigned*, const ArraySpec&, const eckit::Parametrisation&);

template class ArrayT<int>;
template class ArrayT<long>;
template class ArrayT<float>;
template class ArrayT<double>;
template class ArrayT<unsigned long>;

}  // namespace array
}  // namespace atlas
