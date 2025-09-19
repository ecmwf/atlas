/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "SparseMatrixStorage.h"

namespace atlas::linalg {


SparseMatrixStorage::SparseMatrixStorage(SparseMatrixStorage&& other) {
    nnz_     = other.nnz_;
    rows_    = other.rows_;
    cols_    = other.cols_;
    outer_   = std::move(other.outer_);
    inner_   = std::move(other.inner_);
    value_   = std::move(other.value_);
    storage_ = std::move(other.storage_);

    // Invalidate other, just in case it is being used (IT SHOULDNT)
    other.nnz_  = 0;
    other.rows_ = 0;
    other.cols_ = 0;
}

SparseMatrixStorage::SparseMatrixStorage(const SparseMatrixStorage& other) {
    nnz_   = other.nnz_;
    rows_  = other.rows_;
    cols_  = other.cols_;
    {
        array::label label{"sparse_matrix.outer"};
        outer_.reset(atlas::array::Array::create(other.outer_->datatype(), atlas::array::make_shape(other.outer_->size())));
    }
    {
        array::label label{"sparse_matrix.inner"};
        inner_.reset(atlas::array::Array::create(other.inner_->datatype(), atlas::array::make_shape(other.inner_->size())));
    }
    {
        array::label label{"sparse_matrix.value"};
        value_.reset(atlas::array::Array::create(other.value_->datatype(), atlas::array::make_shape(other.value_->size())));
    }
    outer_->copy(*other.outer_);
    inner_->copy(*other.inner_);
    value_->copy(*other.value_);
}

SparseMatrixStorage& SparseMatrixStorage::operator=(SparseMatrixStorage&& other) {
    nnz_     = other.nnz_;
    rows_    = other.rows_;
    cols_    = other.cols_;
    outer_   = std::move(other.outer_);
    inner_   = std::move(other.inner_);
    value_   = std::move(other.value_);
    storage_ = std::move(other.storage_);

    // Invalidate other, just in case it is being used (IT SHOULDNT)
    other.nnz_  = 0;
    other.rows_ = 0;
    other.cols_ = 0;
    return *this;
}


SparseMatrixStorage& SparseMatrixStorage::operator=(const SparseMatrixStorage& other) {
    nnz_   = other.nnz_;
    rows_  = other.rows_;
    cols_  = other.cols_;
    {
        array::label label{"sparse_matrix.outer"};
        outer_.reset(atlas::array::Array::create(other.outer_->datatype(), atlas::array::make_shape(other.outer_->size())));
    }
    {
        array::label label{"sparse_matrix.inner"};
        inner_.reset(atlas::array::Array::create(other.inner_->datatype(), atlas::array::make_shape(other.inner_->size())));
    }
    {
        array::label label{"sparse_matrix.value"};
        value_.reset(atlas::array::Array::create(other.value_->datatype(), atlas::array::make_shape(other.value_->size())));
    }
    outer_->copy(*other.outer_);
    inner_->copy(*other.inner_);
    value_->copy(*other.value_);
    return *this;
}

void SparseMatrixStorage::clear() {
    nnz_ = 0;
    rows_ = 0;
    cols_ = 0;
    storage_.reset();
    outer_.reset();
    inner_.reset();
    value_.reset();
}



void SparseMatrixStorage::updateDevice() const {
    outer_->updateDevice();
    inner_->updateDevice();
    value_->updateDevice();
}

void SparseMatrixStorage::updateHost() const {
    outer_->updateHost();
    inner_->updateHost();
    value_->updateHost();
}

void SparseMatrixStorage::syncDevice() const {
    outer_->syncDevice();
    inner_->syncDevice();
    value_->syncDevice();
}

void SparseMatrixStorage::syncHost() const {
    outer_->syncHost();
    inner_->syncHost();
    value_->syncHost();
}

bool SparseMatrixStorage::hostNeedsUpdate() const {
    return outer_->hostNeedsUpdate() ||
           inner_->hostNeedsUpdate() ||
           value_->hostNeedsUpdate();
}

bool SparseMatrixStorage::deviceNeedsUpdate() const {
    return outer_->deviceNeedsUpdate() ||
           inner_->deviceNeedsUpdate() ||
           value_->deviceNeedsUpdate();
}

void SparseMatrixStorage::setHostNeedsUpdate(bool v) const {
    outer_->setHostNeedsUpdate(v);
    inner_->setHostNeedsUpdate(v);
    value_->setHostNeedsUpdate(v);
}

void SparseMatrixStorage::setDeviceNeedsUpdate(bool v) const {
    outer_->setDeviceNeedsUpdate(v);
    inner_->setDeviceNeedsUpdate(v);
    value_->setDeviceNeedsUpdate(v);
}

bool SparseMatrixStorage::deviceAllocated() const {
    return outer_->deviceAllocated() &&
           inner_->deviceAllocated() &&
           value_->deviceAllocated();
}

void SparseMatrixStorage::allocateDevice() const {
    outer_->allocateDevice();
    inner_->allocateDevice();
    value_->allocateDevice();
}

void SparseMatrixStorage::deallocateDevice() const {
    outer_->deallocateDevice();
    inner_->deallocateDevice();
    value_->deallocateDevice();
}

namespace detail {
    template<typename OutputT, typename InputT>
    void host_copy(const InputT* input_data, array::Array& output) {
        auto size = output.size();
        OutputT* output_data = output.host_data<OutputT>();
        std::copy( input_data, input_data + size, output_data );
    }

    template<typename InputT, typename OutputT>
    void host_copy(const array::Array& input, array::Array& output) {
        host_copy<OutputT>( input.host_data<InputT>(), output );
    }

    template<typename OutputT>
    void host_copy(const array::Array& input, array::Array& output) {
        switch(input.datatype().kind()) {
            case DataType::kind<int>():           return host_copy<int,OutputT>( input, output );
            case DataType::kind<long>():          return host_copy<long,OutputT>( input, output );
            case DataType::kind<float>():         return host_copy<float,OutputT>( input, output );
            case DataType::kind<double>():        return host_copy<double,OutputT>( input, output );
            case DataType::kind<unsigned int>():  return host_copy<unsigned int,OutputT>( input, output );
            case DataType::kind<unsigned long>(): return host_copy<unsigned long,OutputT>( input, output );
            default:  ATLAS_NOTIMPLEMENTED;
        }
    }
    void host_copy(const array::Array& input, array::Array& output) {
        switch(output.datatype().kind()) {
            case DataType::kind<int>():           return host_copy<int>( input, output );
            case DataType::kind<long>():          return host_copy<long>( input, output );
            case DataType::kind<float>():         return host_copy<float>( input, output );
            case DataType::kind<double>():        return host_copy<double>( input, output );
            case DataType::kind<unsigned int>():  return host_copy<unsigned int>( input, output );
            case DataType::kind<unsigned long>(): return host_copy<unsigned long>( input, output );
            default:  ATLAS_NOTIMPLEMENTED;
        }
    }
}



}