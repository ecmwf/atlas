/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file Table.h
/// @author Willem Deconinck
/// @date October 2015
/// This file contains the classes
///   - Table
///
/// It is important to note that connectivity access functions are
/// inlined for optimal performance. The connectivity itself is internally
/// stored with base 1, to be compatible with Fortran access.
/// C++ access operators however convert the resulting connectivity to base 0.

#pragma once

#include <type_traits>

#include "atlas/array/Table.h"
#include "atlas/library/config.h"

namespace atlas {
namespace array {

// --------------------------------------------------------------------------

#if ATLAS_HAVE_FORTRAN
#define INDEX_REF Index
#define FROM_FORTRAN -1
#define TO_FORTRAN +1
#else
#define INDEX_REF *
#define FROM_FORTRAN
#define TO_FORTRAN
#endif

namespace detail {

// TableIndex:
// Helper class that does +1 and -1 operations on stored values

class TableIndex {
public:
    enum
    {
        BASE = 1
    };

public:
    TableIndex(idx_t* idx): idx_(idx) {}
    void set(const idx_t& value) { *(idx_) = value + BASE; }
    idx_t get() const { return *(idx_)-BASE; }
    void operator       =(const idx_t& value) { set(value); }
    TableIndex& operator=(const TableIndex& other) {
        set(other.get());
        return *this;
    }
    TableIndex& operator+(const idx_t& value) {
        *(idx_) += value;
        return *this;
    }
    TableIndex& operator-(const idx_t& value) {
        *(idx_) -= value;
        return *this;
    }
    TableIndex& operator--() {
        --(*(idx_));
        return *this;
    }
    TableIndex& operator++() {
        ++(*(idx_));
        return *this;
    }
    TableIndex& operator+=(const idx_t& value) {
        *(idx_) += value;
        return *this;
    }
    TableIndex& operator-=(const idx_t& value) {
        *(idx_) -= value;
        return *this;
    }

    // implicit conversion
    operator idx_t() const { return get(); }

private:
    idx_t* idx_;
};

}  // namespace detail

// ------------------------------------------------------------------------------------------------------

template <bool ReadOnly>
class TableRow {
#if ATLAS_HAVE_FORTRAN
    typedef detail::TableIndex Index;
#else
    typedef idx_t Index;
#endif

    using ReturnType = typename std::conditional<ReadOnly, idx_t, Index>::type;

public:
    ATLAS_HOST_DEVICE
    TableRow(const idx_t* data, size_t size): data_(const_cast<idx_t*>(data)), size_(size) {}

    ATLAS_HOST_DEVICE
    idx_t operator()(size_t i) const { return data_[i] FROM_FORTRAN; }

    ATLAS_HOST_DEVICE
    ReturnType operator()(size_t i) { return INDEX_REF(data_ + i); }

    ATLAS_HOST_DEVICE
    size_t size() const { return size_; }

    // TODO: Following function should only be allowed to compile if
    // ReadOnly=false (SFINAE?)
    TableRow& operator=(const idx_t column_values[]) {
        assert(not ReadOnly);
        for (size_t n = 0; n < size_; ++n) {
            data_[n] = column_values[n] TO_FORTRAN;
        }
        return *this;
    }

private:
    idx_t* data_;
    size_t size_;
};

// ------------------------------------------------------------------------------------------------------

template <bool ReadOnly>
class TableView : public util::Object {
#if ATLAS_HAVE_FORTRAN
    using Index = typename std::conditional<ReadOnly, idx_t, detail::TableIndex>::type;
#else
    using Index = idx_t;
#endif

public:
    typedef TableRow<ReadOnly> Row;
    typedef TableRow<true> ConstRow;

    static constexpr unsigned short _values_ = 0;
    static constexpr unsigned short _displs_ = 1;
    static constexpr unsigned short _counts_ = 2;

public:
    using value_type = idx_t;
    using Data       = typename std::conditional<ReadOnly, const idx_t, idx_t>::type;

    template <typename ReturnType, bool ReadOnly_>
    struct Access_t {};

    template <typename ReturnType>
    struct Access_t<ReturnType, true> {
        Access_t(const TableView<ReadOnly>* tv): tv_(const_cast<TableView<ReadOnly>*>(tv)) {}
        TableView<ReadOnly>* tv_;
        idx_t apply(size_t row, size_t col) const { return tv_->values_(tv_->displs_(row) + col) FROM_FORTRAN; }
    };

    template <typename ReturnType>
    struct Access_t<ReturnType, false> {
        Access_t(const TableView<ReadOnly>* tv): tv_(const_cast<TableView<ReadOnly>*>(tv)) {}
        TableView<ReadOnly>* tv_;
        Index apply(size_t row, size_t col) const { return Index(&tv_->values_(tv_->displs_(row) + col)); }
    };

    using Access      = Access_t<Index, ReadOnly>;
    using ConstAccess = Access_t<idx_t, true>;

    //-- Constructors

    TableView(const Table& table, bool host = true);

    TableView(const TableView& other);

    TableView operator=(const TableView& other);

    ~TableView() {}

    //-- Accessors

    /// @brief Number of rows in the connectivity table
    ATLAS_HOST_DEVICE
    size_t rows() const { return rows_; }

    /// @brief Number of columns for specified row in the connectivity table
    ATLAS_HOST_DEVICE
    size_t cols(size_t row_idx) const { return counts_(row_idx); }

    /// @brief Maximum value for number of columns over all rows
    ATLAS_HOST_DEVICE
    size_t maxcols() const { return maxcols_; }

    /// @brief Minimum value for number of columns over all rows
    ATLAS_HOST_DEVICE
    size_t mincols() const { return mincols_; }

    /// @brief Access to connectivity table elements for given row and column
    /// The returned index has base 0 regardless if ATLAS_HAVE_FORTRAN is defined.
    ATLAS_HOST_DEVICE
    idx_t operator()(size_t row_idx, size_t col_idx) const {
        assert(counts_(row_idx) > col_idx);
        return const_access_.apply(row_idx, col_idx);
    }

    /// @brief Access to connectivity table elements for given row and column
    /// The returned index has base 0 regardless if ATLAS_HAVE_FORTRAN is defined.
    ATLAS_HOST_DEVICE
    Index operator()(size_t row_idx, size_t col_idx) {
        assert(counts_(row_idx) > col_idx);
        return access_.apply(row_idx, col_idx);
    }

    /// @brief Access to raw data.
    /// Note that the connectivity base is 1 in case ATLAS_HAVE_FORTRAN is
    /// defined.
    const idx_t* data() const { return values_.data(); }
    Data* data() { return values_.data(); }

    ATLAS_HOST_DEVICE
    size_t size() const { return values_.size(); }

    ATLAS_HOST_DEVICE
    idx_t missing_value() const { return missing_value_; }

    ATLAS_HOST_DEVICE
    ConstRow row(size_t row_idx) const { return ConstRow(values_.data() + displs_(row_idx), counts_(row_idx)); }

    ATLAS_HOST_DEVICE
    Row row(size_t row_idx) { return Row(values_.data() + displs_(row_idx), counts_(row_idx)); }

    ///-- Modifiers

    ATLAS_HOST_DEVICE
    size_t displs(const size_t row) const { return displs_(row); }

private:
    const size_t* displs() const { return displs_.data(); }
    const size_t* counts() const { return counts_.data(); }

private:
    bool host_;
    idx_t missing_value_;
    size_t rows_;
    size_t maxcols_;
    size_t mincols_;
    ArrayView<idx_t, 1> values_;
    ArrayView<size_t, 1> displs_;
    ArrayView<size_t, 1> counts_;
    ConstAccess const_access_;
    Access access_;
};

// -----------------------------------------------------------------------------------------------------

#undef FROM_FORTRAN
#undef TO_FORTRAN
#undef INDEX_REF

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
