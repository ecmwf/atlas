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

#pragma once

#include <array>
#include <type_traits>

#include "atlas/util/Object.h"

#include "atlas/array.h"
#include "atlas/library/config.h"

namespace atlas {
namespace array {

// --------------------------------------------------------------------------

/// @brief Table
/// @author Willem Deconinck
///
/// Container for tables. This is e.g.
/// for a node-to-X connectivity.
///   connectivity = [
///     1   2   3   4   5   6       # node 1
///     7   8                       # node 2
///     9  10  11  12               # node 3
///    13  14  15                   # node 4
///    16  17  18  19  20  21       # node 5
///   ]
/// There are 2 modes of construction:
///   - It wraps existing raw data
///   - It owns the connectivity data
///
/// In case ATLAS_HAVE_FORTRAN is defined (which is usually the case),
/// the raw data will be stored with base 1 for Fortran interoperability.
/// The operator(row,col) in C++ will then do the conversion to base 0.
///
/// In the first mode of construction, the connectivity table cannot be resized.
/// In the second mode of construction, resizing is possible

class Table : public util::Object {
private:
    static constexpr unsigned short _values_ = 0;
    static constexpr unsigned short _displs_ = 1;
    static constexpr unsigned short _counts_ = 2;

public:
    //-- Constructors

    /// @brief Construct connectivity table that needs resizing a-posteriori
    /// Data is owned
    Table(const std::string& name = "");

    ~Table();

private:
    /// @brief Construct connectivity table wrapping existing raw data.
    /// No resizing can be performed as data is not owned.
    Table(idx_t values[], size_t rows, size_t displs[], size_t counts[]);

public:
    //-- Accessors

    /// @brief Name associated to this Connetivity
    const std::string& name() const { return name_; }

    /// @brief Rename this Connectivity
    void rename(const std::string& name) { name_ = name; }

    /// @brief Number of rows in the connectivity table
    size_t rows() const { return rows_; }

    /// @brief Number of columns for specified row in the connectivity table
    size_t cols(size_t row_idx) const { return counts_(row_idx); }

    /// @brief Maximum value for number of columns over all rows
    size_t maxcols() const { return maxcols_; }

    /// @brief Minimum value for number of columns over all rows
    size_t mincols() const { return mincols_; }

    /// @brief Value that is given to unassigned entries
    idx_t missing_value() const { return missing_value_; }

    /// @brief Number of values stored in the table
    size_t size() const { return data_[_values_]->size(); }

    /// @brief Return memory footprint of table
    virtual size_t footprint() const;

    /// @brief Clone data to device
    virtual void updateDevice() const;

    /// @brief Clone data from device
    virtual void updateHost() const;

    /// @brief Synchronise data between host and device
    virtual void syncHostDevice() const;

    /// @brief Check if data is valid
    virtual bool valid() const;

    /// @brief Check if data is present on host
    virtual bool hostNeedsUpdate() const;

    /// @brief Check if data is present on device
    virtual bool deviceNeedsUpdate() const;

    /// @brief Print all values unformatted to output stream
    void dump(std::ostream&) const;

    /// @brief Check if data is owned or wrapped
    bool owns() { return owns_; }

    ///-- Modifiers

    /// @brief Resize connectivity, and add given rows
    /// @note Can only be used when data is owned.
    virtual void add(size_t rows, size_t cols, const idx_t values[], bool fortran_array = false);

    /// @brief Resize connectivity, and add given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void add(size_t rows, size_t cols);

    /// @brief Resize connectivity, and add given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void add(size_t rows, const size_t cols[]);

    /// @brief Resize connectivity, and insert given rows
    /// @note Can only be used when data is owned.
    virtual void insert(size_t position, size_t rows, size_t cols, const idx_t values[], bool fortran_array = false);

    /// @brief Resize connectivity, and insert given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void insert(size_t position, size_t rows, size_t cols);

    /// @brief Resize connectivity, and insert given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void insert(size_t position, size_t rows, const size_t cols[]);

    /// @brief Resize connectivity, and insert given rows with missing values
    /// @note Invalidates non-owned Table
    virtual void clear();

private:
    ///-- Internal helper functions

    void resize_values(size_t old_size, size_t size, bool initialize, const idx_t values[], bool fortran_array);

    void resize_counts_and_displs(size_t size);

    void insert_counts_and_displs(size_t position, size_t rows);

    void resize_values(size_t size);

    void insert_values(size_t position, size_t size);

private:
    template <bool ReadOnly>
    friend class TableView;

    std::string name_;
    bool owns_;
    std::array<array::Array*, 3> data_;
    idx_t missing_value_;
    size_t rows_;
    size_t maxcols_;
    size_t mincols_;
    ArrayView<size_t, 1> displs_;
    ArrayView<size_t, 1> counts_;
    ArrayView<idx_t, 1> values_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
