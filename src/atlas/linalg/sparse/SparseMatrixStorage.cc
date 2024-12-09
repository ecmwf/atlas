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
    outer_.reset(atlas::array::Array::create(other.outer_->datatype(), atlas::array::make_shape(other.outer_->size())));
    inner_.reset(atlas::array::Array::create(other.inner_->datatype(), atlas::array::make_shape(other.inner_->size())));
    value_.reset(atlas::array::Array::create(other.value_->datatype(), atlas::array::make_shape(other.value_->size())));
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
    outer_.reset(atlas::array::Array::create(other.outer_->datatype(), atlas::array::make_shape(other.outer_->size())));
    inner_.reset(atlas::array::Array::create(other.inner_->datatype(), atlas::array::make_shape(other.inner_->size())));
    value_.reset(atlas::array::Array::create(other.value_->datatype(), atlas::array::make_shape(other.value_->size())));
    outer_->copy(*other.outer_);
    inner_->copy(*other.inner_);
    value_->copy(*other.value_);
    return *this;
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




}