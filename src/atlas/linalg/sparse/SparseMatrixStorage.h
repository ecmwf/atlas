/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <algorithm>
#include <any>
#include <memory>
#include <utility>

#include "eckit/linalg/SparseMatrix.h"

#include "atlas/array.h"
#include "atlas/runtime/Exception.h"
#include "atlas/linalg/sparse/SparseMatrixView.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

/// @brief SparseMatrixStorage
/// Storage class based on atlas::array::Array with GPU offload capability
/// This class contains only the data and no operators, iterators, or special
/// sparse matrix construction from triplets etc.
///
/// The construction is handled via SparseMatrixConvertor class.
/// Examples:
///
/// To construct it, taking ownership of a constructed Eigen::SparseMatrix
///
///     SparseMatrixStorage s = make_sparse_matrix_storage(std::move(eigen_matrix));
///
/// To construct it, taking a copy of a constructed Eigen::SparseMatrix
///
///     SparseMatrixStorage s = make_sparse_matrix_storage(eigen_matrix);
///
/// To construct it, taking ownership of a constructed eckit::linalg::SparseMatrix, avoiding copies if data types match
///
///     SparseMatrixStorage s = make_sparse_matrix_storage(std::move(eckit_matrix));
///
///
/// To construct it, taking a copy of a constructed eckit::linalg::SparseMatrix (no std::move)
///
///     SparseMatrixStorage s = make_sparse_matrix_storage(eckit_matrix);
///
///
/// It is also possible to initialise empty and move into it at later stage:
///
///     SparseMatrixStorage s;
///     s = make_sparse_matrix_storage(eckit_matrix);
///
///
/// To construct it, taking a single precision copy of a constructed eckit::linalg::SparseMatrix
///
///     s = make_sparse_matrix_storage<float>(eckit_matrix);
///

class SparseMatrixStorage {
public:

    // Virtual destructor
    virtual ~SparseMatrixStorage() = default;

    /// Default constructor
    SparseMatrixStorage() = default;

    /// Move constructor, takes ownership!
    SparseMatrixStorage(SparseMatrixStorage&& other);

    /// Copy constructor, makes copy!
    SparseMatrixStorage(const SparseMatrixStorage& other);

    /// Copy from a SparseMatrixView
    template<typename Value, typename Index>
    SparseMatrixStorage(const SparseMatrixView<Value,Index>& host_view);

    /// Move assign from other SparseMatrixStorage, takes ownership!
    SparseMatrixStorage& operator=(SparseMatrixStorage&& other);

    /// Copy assign from other SparseMatrixStorage, takes ownership!
    SparseMatrixStorage& operator=(const SparseMatrixStorage& other);

    /// Empty if rows and cols are zero.
    bool empty() const { return rows_ == 0 && cols_ == 0; }

    /// Footprint in bytes in host memory space
    std::size_t footprint() const { return value_->footprint() + outer_->footprint() + inner_->footprint(); }

    const atlas::array::Array& value() const { return *value_; }
    const atlas::array::Array& outer() const { return *outer_; }
    const atlas::array::Array& inner() const { return *inner_; }
    std::size_t rows() const { return rows_;}
    std::size_t cols() const { return cols_;}
    std::size_t nnz()  const { return nnz_;}

    void updateDevice() const;

    void updateHost() const;

    bool hostNeedsUpdate() const;

    bool deviceNeedsUpdate() const;

    void setHostNeedsUpdate(bool v) const;

    void setDeviceNeedsUpdate(bool v) const;

    bool deviceAllocated() const;

    void allocateDevice() const;

    void deallocateDevice() const;

    static SparseMatrixStorage make(
        std::size_t rows,
        std::size_t cols,
        std::size_t nnz,
        std::unique_ptr<atlas::array::Array>&& value,
        std::unique_ptr<atlas::array::Array>&& inner,
        std::unique_ptr<atlas::array::Array>&& outer,
        std::any&& storage) {
            SparseMatrixStorage S;
            S.rows_    = rows;
            S.cols_    = cols;
            S.nnz_     = nnz;
            S.outer_   = std::move(outer);
            S.inner_   = std::move(inner);
            S.value_   = std::move(value);
            S.storage_ = std::move(storage);
            return S;
    }

    void swap(SparseMatrixStorage& other) {
        std::swap(other.rows_,    rows_);
        std::swap(other.cols_,    cols_);
        std::swap(other.nnz_,     nnz_);
        std::swap(other.value_,   value_);
        std::swap(other.inner_,   inner_);
        std::swap(other.outer_,   outer_);
        std::swap(other.storage_, storage_);
    }

    bool contains(DataType value_type, DataType index_type) {
        if (empty()) {
            return false;
        }
        return value_->datatype() == value_type && outer_->datatype() == index_type;
    }

    template <typename value_type, typename index_type>
    bool contains() {
        return contains(make_datatype<value_type>(), make_datatype<index_type>());
    }
    
protected:
    std::size_t nnz_{0};
    std::size_t rows_{0};
    std::size_t cols_{0};
    std::unique_ptr<atlas::array::Array> outer_;
    std::unique_ptr<atlas::array::Array> inner_;
    std::unique_ptr<atlas::array::Array> value_;

    std::any storage_;
        // This storage is usually empty.
        // It is used to allow to move alternative sparse matrix formats into it,
        // and then wrap this data using the atlas Arrays

private:
    template<typename ValueT>
    void host_copy(const ValueT* input_data, atlas::array::Array& output) {
        auto size = output.size();
        ValueT* output_data = output.host_data<ValueT>();
        std::copy( input_data, input_data + size, output_data );
    }
};

//----------------------------------------------------------------------------------------------------------------------

template<typename Value, typename Index>
SparseMatrixStorage::SparseMatrixStorage(const SparseMatrixView<Value,Index>& host_view) {
    nnz_   = host_view.nnz();
    rows_  = host_view.rows();
    cols_  = host_view.cols();
    outer_.reset(atlas::array::Array::create<Index>(host_view.outer_size()));
    inner_.reset(atlas::array::Array::create<Index>(host_view.inner_size()));
    value_.reset(atlas::array::Array::create<Value>(host_view.value_size()));
    host_copy(host_view.outer(),*outer_);
    host_copy(host_view.inner(),*inner_);
    host_copy(host_view.value(),*value_);
}

//----------------------------------------------------------------------------------------------------------------------

namespace detail {
    void host_copy(const array::Array& input, array::Array& output);
}

template<typename value_type, typename index_type = eckit::linalg::Index>
SparseMatrixStorage make_sparse_matrix_storage(const SparseMatrixStorage& other) {
    auto rows  = other.rows();
    auto cols  = other.cols();
    auto nnz   = other.nnz();
    std::unique_ptr<array::Array> value(array::Array::create<value_type>(nnz));
    std::unique_ptr<array::Array> inner(array::Array::create<index_type>(nnz));
    std::unique_ptr<array::Array> outer(array::Array::create<index_type>(rows+1));
    detail::host_copy(other.value(), *value);
    detail::host_copy(other.inner(), *inner);
    detail::host_copy(other.outer(), *outer);
    return SparseMatrixStorage::make(rows,cols,nnz,std::move(value), std::move(inner), std::move(outer), std::any());
}

template<typename value_type, typename index_type = eckit::linalg::Index>
SparseMatrixStorage make_sparse_matrix_storage(SparseMatrixStorage&& other) {
    SparseMatrixStorage S;

    if (other.contains<value_type,index_type>()) {
        S = std::move(other);
    }
    else {
        auto rows  = other.rows();
        auto cols  = other.cols();
        auto nnz   = other.nnz();
        std::unique_ptr<array::Array> value(array::Array::create<value_type>(nnz));
        std::unique_ptr<array::Array> inner(array::Array::create<index_type>(nnz));
        std::unique_ptr<array::Array> outer(array::Array::create<index_type>(rows+1));
        detail::host_copy(other.value(), *value);
        detail::host_copy(other.inner(), *inner);
        detail::host_copy(other.outer(), *outer);
        S = SparseMatrixStorage::make(rows,cols,nnz,std::move(value), std::move(inner), std::move(outer), std::any());
    }
    return S;
}

//----------------------------------------------------------------------------------------------------------------------

template<typename Value, typename Index = eckit::linalg::Index>
inline SparseMatrixView<Value,Index> make_host_view(const SparseMatrixStorage& m) {
    if(m.rows() == 0 && m.cols() == 0) {
        return SparseMatrixView<Value,Index>();
    }
    if( m.value().datatype().kind() != DataType::kind<Value>() || m.outer().datatype().kind() != DataType::kind<Index>() ) {
        ATLAS_THROW_EXCEPTION("Cannot make_host_view<" + DataType::str<Value>() + "," << DataType::str<Index>() +
            ">(const SparseMatrixStorage&) from SparseMatrixStorage containing values of type <" + m.value().datatype().str() + "> and indices of type <" + m.outer().datatype().str() +">" );
    }
    return SparseMatrixView<Value,Index> {
        m.rows(),
        m.cols(),
        m.nnz(),
        m.value().host_data<Value>(),
        m.inner().host_data<Index>(),
        m.outer().host_data<Index>()
    };
}

//----------------------------------------------------------------------------------------------------------------------

template<typename Value, typename Index = eckit::linalg::Index>
inline SparseMatrixView<Value,Index> make_device_view(const SparseMatrixStorage& m) {
    if(m.rows() == 0 && m.cols() == 0) {
        return SparseMatrixView<Value,Index>();
    }
    if( m.value().datatype().kind() != DataType::kind<Value>() || m.outer().datatype().kind() != DataType::kind<Index>() ) {
        ATLAS_THROW_EXCEPTION("Cannot make_device_view<" + DataType::str<Value>() + "," << DataType::str<Index>() +
            ">(const SparseMatrixStorage&) from SparseMatrixStorage containing values of type <" + m.value().datatype().str() + "> and indices of type <" + m.outer().datatype().str() +">" );
    }
    if( ! m.deviceAllocated() ) {
        ATLAS_THROW_EXCEPTION("Cannot make_device_view(const SparseMatrixStorage&) as the device data is not allocated");
    }
    return SparseMatrixView<Value,Index>{
        m.rows(),
        m.cols(),
        m.nnz(),
        m.value().device_data<Value>(),
        m.inner().device_data<Index>(),
        m.outer().device_data<Index>()
    };
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace linalg
}  // namespace atlas