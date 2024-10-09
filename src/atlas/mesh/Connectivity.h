/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file Connectivity.h
/// @author Willem Deconinck
/// @date October 2015
/// @details
/// This file contains Connectivity classes
///   - IrregularConnectivity
///   - BlockConnectivity
///   - MultiBlockConnectivity
///
/// It is important to note that connectivity access functions are
/// inlined for optimal performance. The connectivity itself is internally
/// stored with base 1, to be compatible with Fortran access.
/// C++ access operators however convert the resulting connectivity to base 0.

#pragma once

#include <array>
#include <cstring>
#include <initializer_list>

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/DataType.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/SVector.h"
#include "atlas/array_fwd.h"
#include "atlas/library/config.h"
#include "atlas/util/Object.h"

namespace eckit {
class Stream;
}

namespace atlas {
namespace mesh {

constexpr size_t MAX_STRING_SIZE() {
    return 60;
}

template <typename ConnectivityImpl>
class ConnectivityInterface : public util::Object, DOXYGEN_HIDE(public ConnectivityImpl) {
    using ConnectivityImpl::ConnectivityImpl;
    using util::Object::Object;
};

// Classes defined in this file:
class IrregularConnectivityImpl;
class BlockConnectivityImpl;
class MultiBlockConnectivityImpl;

// --------------------------------------------------------------------------

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if ATLAS_HAVE_FORTRAN
#define INDEX_DEREF Index
#define INDEX_REF Index
#define FROM_FORTRAN -1
#define TO_FORTRAN +1
#else
#define INDEX_DEREF *
#define INDEX_REF idx_t&
#define FROM_FORTRAN
#define TO_FORTRAN
#endif
#endif

/// @brief Irregular Connectivity Table
/// @author Willem Deconinck
///
/// Container for irregular connectivity tables. This is e.g. the case
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
/// The operator(row,col) will then do the conversion to base 0.
///
/// In the first mode of construction, the connectivity table cannot be resized.
/// In the second mode of construction, resizing is possible

namespace detail {
// FortranIndex:
// Helper class that does +1 and -1 operations on stored values

class ConnectivityIndex {
public:
    enum
    {
        BASE = 1
    };

public:
    ATLAS_HOST_DEVICE ConnectivityIndex(const ConnectivityIndex& other) { set(other.get()); }
    ATLAS_HOST_DEVICE ConnectivityIndex(idx_t* idx): idx_(idx) {}
    ATLAS_HOST_DEVICE void set(const idx_t& value) { *(idx_) = value + BASE; }
    ATLAS_HOST_DEVICE idx_t get() const { return *(idx_)-BASE; }
    ATLAS_HOST_DEVICE void operator=(const idx_t& value) { set(value); }
    ATLAS_HOST_DEVICE ConnectivityIndex& operator+(const idx_t& value) {
        *(idx_) += value;
        return *this;
    }
    ATLAS_HOST_DEVICE ConnectivityIndex& operator-(const idx_t& value) {
        *(idx_) -= value;
        return *this;
    }
    ATLAS_HOST_DEVICE ConnectivityIndex& operator--() {
        --(*(idx_));
        return *this;
    }
    ATLAS_HOST_DEVICE ConnectivityIndex& operator++() {
        ++(*(idx_));
        return *this;
    }

    // implicit conversion
    ATLAS_HOST_DEVICE operator idx_t() const { return get(); }

private:
    idx_t* idx_;
};
}  // namespace detail

class ConnectivityRow {
#if ATLAS_HAVE_FORTRAN
    using Index = detail::ConnectivityIndex;
#else
    using Index = idx_t&;
#endif
public:
    ATLAS_HOST_DEVICE
    ConnectivityRow(idx_t* data, idx_t size): data_(data), size_(size) {}

    template <typename Int>
    ATLAS_HOST_DEVICE idx_t operator()(Int i) const {
        return data_[i] FROM_FORTRAN;
    }

    template <typename Int>
    ATLAS_HOST_DEVICE INDEX_REF operator()(Int i) {
        return INDEX_DEREF(data_ + i);
    }

    ATLAS_HOST_DEVICE
    idx_t size() const { return size_; }

private:
    idx_t* data_;
    idx_t size_;
};

class IrregularConnectivityImpl {
public:
    typedef ConnectivityRow Row;

public:
    //-- Constructors

    /// @brief Construct connectivity table that needs resizing a-posteriori
    /// Data is owned
    IrregularConnectivityImpl(const std::string& name = "");

    /// @brief Construct connectivity table wrapping existing raw data.
    /// No resizing can be performed as data is not owned.
    IrregularConnectivityImpl(idx_t values[], idx_t rows, idx_t displs[], idx_t counts[]);

#if ATLAS_HIC_COMPILER
    /// @brief Copy ctr (only to be used when calling a cuda kernel)
    // This ctr has to be defined in the header, since __CUDACC__ will identify
    // whether
    // it is compiled it for a GPU kernel
    IrregularConnectivityImpl(const IrregularConnectivityImpl& other):
        owns_(false),
        values_(other.values_),
        displs_(other.displs_),
        counts_(other.counts_),
        missing_value_(other.missing_value_),
        rows_(other.rows_),
        maxcols_(other.maxcols_),
        mincols_(other.mincols_),
        ctxt_(nullptr) {}
#endif

    /// @brief Construct a mesh from a Stream (serialization)
    explicit IrregularConnectivityImpl(eckit::Stream&);

    virtual ~IrregularConnectivityImpl();

    //-- Accessors

    /// @brief Name associated to this Connetivity
    const std::string name() const { return std::string(name_); }

    /// @brief Rename this Connectivity
    void rename(const std::string& name) {
        std::strncpy(name_, name.c_str(), std::max(name.size(), MAX_STRING_SIZE()));
    }

    /// @brief Number of rows in the connectivity table
    ATLAS_HOST_DEVICE
    idx_t rows() const { return rows_; }

    /// @brief Number of columns for specified row in the connectivity table
    ATLAS_HOST_DEVICE
    idx_t cols(idx_t row_idx) const { return counts_[row_idx]; }

    /// @brief Maximum value for number of columns over all rows
    ATLAS_HOST_DEVICE
    idx_t maxcols() const { return maxcols_; }

    /// @brief Minimum value for number of columns over all rows
    ATLAS_HOST_DEVICE
    idx_t mincols() const { return mincols_; }

    /// @brief Access to connectivity table elements for given row and column
    /// The returned index has base 0 regardless if ATLAS_HAVE_FORTRAN is defined.
    ATLAS_HOST_DEVICE
    idx_t operator()(idx_t row_idx, idx_t col_idx) const;

    idx_t size() const { return values_.size(); }

    ATLAS_HOST_DEVICE
    idx_t missing_value() const { return missing_value_; }

    ATLAS_HOST_DEVICE
    Row row(idx_t row_idx) const;

    // -- Modifiers

    /// @brief Modify row with given values. Values must be given with base 0
    void set(idx_t row_idx, const idx_t column_values[]);

    /// @brief Modify (row,col) with given value. Value must be given with base 0
    void set(idx_t row_idx, idx_t col_idx, const idx_t value);

    /// @brief Resize connectivity
    /// @note Can only be used when data is owned.
    virtual void resize(idx_t old_size, idx_t size, bool initialize, const idx_t values[], bool fortran_array);

    /// @brief Resize connectivity, and add given rows
    /// @note Can only be used when data is owned.
    virtual void add(idx_t rows, idx_t cols, const idx_t values[], bool fortran_array = false);

    /// @brief Resize connectivity, and add given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void add(idx_t rows, idx_t cols);

    /// @brief Resize connectivity, and add given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void add(idx_t rows, const idx_t cols[]);

    /// @brief Resize connectivity, and copy from a BlockConnectivity
    virtual void add(const BlockConnectivityImpl& block);

    /// @brief Resize connectivity, and insert given rows
    /// @note Can only be used when data is owned.
    virtual void insert(idx_t position, idx_t rows, idx_t cols, const idx_t values[], bool fortran_array = false);

    /// @brief Resize connectivity, and insert given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void insert(idx_t position, idx_t rows, idx_t cols);

    /// @brief Resize connectivity, and insert given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void insert(idx_t position, idx_t rows, const idx_t cols[]);

    virtual void clear();

    virtual size_t footprint() const;

    idx_t displs(const idx_t row) const { return displs_[row]; }

    void dump(std::ostream& os) const;
    friend std::ostream& operator<<(std::ostream& os, const IrregularConnectivityImpl& p) {
        p.dump(os);
        return os;
    }

    virtual void encode(eckit::Stream& s) const { encode_(s); }
    virtual void decode(eckit::Stream& s) { decode_(s); }

protected:
    bool owns() { return owns_; }
    const idx_t* displs() const { return displs_.data(); }
    const idx_t* counts() const { return counts_.data(); }

    /// @brief Serialization to Stream
    void encode_(eckit::Stream&) const;

    /// @brief Serialization from Stream
    void decode_(eckit::Stream&);

    friend eckit::Stream& operator<<(eckit::Stream& s, const IrregularConnectivityImpl& x) {
        x.encode(s);
        return s;
    }

    friend eckit::Stream& operator>>(eckit::Stream& s, IrregularConnectivityImpl& x) {
        x.decode(s);
        return s;
    }

private:
    void on_delete();
    void on_update();

private:
    char name_[MAX_STRING_SIZE()];
    bool owns_;

protected:
    array::SVector<idx_t> values_;
    array::SVector<idx_t> displs_;
    array::SVector<idx_t> counts_;

    idx_t missing_value_;
    idx_t rows_;
    idx_t maxcols_;
    idx_t mincols_;

public:
    typedef void* ctxt_t;
    typedef void (*callback_t)(ctxt_t);

private:
    friend class ConnectivityPrivateAccess;
    ctxt_t ctxt_;
    callback_t callback_update_;
    callback_t callback_delete_;
};

// ----------------------------------------------------------------------------------------------

/// @brief Connectivity contiguously composed of multiple BlockConnectivityImpl
/// @author Willem Deconinck
///
/// Container for connectivity tables that are layed out in memory as multiple
/// BlockConnectivities stitched together.
/// This is e.g. the case for a element-node connectivity, with element-types
/// grouped together:
/// @code{.sh}
///  connectivity = [
///     # triangles  (block 0)
///         1   2   3
///         4   5   6
///         7   8   9
///        10  11  12
///     # quadrilaterals (block 1)
///        13  14  15  16
///        17  18  19  20
///        21  22  23  24
///  ]
/// @endcode
///
/// This class can also be interpreted as a atlas::IrregularConnectivity without
/// distinction between blocks
///
/// There are 2 modes of construction:
///   - It wraps existing raw data
///   - It owns the connectivity data
///
/// In case ATLAS_HAVE_FORTRAN is defined (which is usually the case),
/// the raw data will be stored with base 1 for Fortran interoperability.
/// The operator(row,col) will then do the conversion to base 0.
///
/// In the first mode of construction, the connectivity table cannot be resized.
/// In the second mode of construction, resizing is possible
class MultiBlockConnectivityImpl : public IrregularConnectivityImpl {
public:
    //-- Constructors

    /// @brief Construct connectivity table that needs resizing a-posteriori
    /// Data is owned
    MultiBlockConnectivityImpl(const std::string& name = "");

    /// @brief Construct a mesh from a Stream (serialization)
    explicit MultiBlockConnectivityImpl(eckit::Stream&);

    virtual ~MultiBlockConnectivityImpl();

    //-- Accessors

    /// @brief Number of blocks
    ATLAS_HOST_DEVICE
    idx_t blocks() const { return blocks_; }

    /// @brief Access to a block connectivity
    ATLAS_HOST_DEVICE
    const BlockConnectivityImpl& block(idx_t block_idx) const { return block_[block_idx]; }
    ATLAS_HOST_DEVICE
    BlockConnectivityImpl& block(idx_t block_idx) { return block_[block_idx]; }

    /// @brief Access to connectivity table elements for given row and column
    /// The row_idx counts up from 0, from block 0, as in IrregularConnectivity
    /// The returned index has base 0 regardless if ATLAS_HAVE_FORTRAN is defined.
    ATLAS_HOST_DEVICE
    idx_t operator()(idx_t row_idx, idx_t col_idx) const;

    /// @brief Access to connectivity table elements for given row and column
    /// The block_row_idx counts up from zero for every block_idx.
    /// The returned index has base 0 regardless if ATLAS_HAVE_FORTRAN is defined.
    ATLAS_HOST_DEVICE
    idx_t operator()(idx_t block_idx, idx_t block_row_idx, idx_t block_col_idx) const;

    //-- Modifiers

    /// @brief Resize connectivity, and add given rows as a new block
    /// @note Can only be used when data is owned.
    virtual void add(idx_t rows, idx_t cols, const idx_t values[], bool fortran_array = false);

    /// @brief Resize connectivity, and copy from a BlockConnectivity to a new
    /// block
    /// @note Can only be used when data is owned.
    virtual void add(const BlockConnectivityImpl&);

    /// @brief Resize connectivity, and add given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void add(idx_t rows, idx_t cols);

    /// @brief Resize connectivity, and add given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void add(idx_t rows, const idx_t cols[]);

    /// @brief Resize connectivity, and insert given rows
    /// @note Can only be used when data is owned.
    virtual void insert(idx_t position, idx_t rows, idx_t cols, const idx_t values[], bool fortran_array = false);

    /// @brief Resize connectivity, and insert given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void insert(idx_t position, idx_t rows, idx_t cols);

    /// @brief Resize connectivity, and insert given rows with missing values
    /// @note Can only be used when data is owned.
    virtual void insert(idx_t position, idx_t rows, const idx_t cols[]);

    virtual void clear();

    virtual size_t footprint() const;

    virtual void encode(eckit::Stream& s) const {
        IrregularConnectivityImpl::encode_(s);
        encode_(s);
    }
    virtual void decode(eckit::Stream& s) {
        IrregularConnectivityImpl::decode_(s);
        decode_(s);
    }

protected:
    /// @brief Serialization to Stream
    void encode_(eckit::Stream&) const;

    /// @brief Serialization from Stream
    void decode_(eckit::Stream&);

    friend eckit::Stream& operator<<(eckit::Stream& s, const MultiBlockConnectivityImpl& x) {
        x.encode(s);
        return s;
    }

    friend eckit::Stream& operator>>(eckit::Stream& s, MultiBlockConnectivityImpl& x) {
        x.decode(s);
        return s;
    }

private:
    void rebuild_block_connectivity();

private:
    idx_t blocks_;
    array::SVector<idx_t> block_displs_;
    array::SVector<idx_t> block_cols_;

    array::SVector<BlockConnectivityImpl> block_;
};

// -----------------------------------------------------------------------------------------------------

/// @brief Block Connectivity table
/// @author Willem Deconinck
/// Container for connectivity tables that are layed out in memory as a block.
/// Every row has the same
/// number of columns.
///
/// There are 2 modes of construction:
///   - It wraps existing raw data
///   - It owns the connectivity data
///
/// In case ATLAS_HAVE_FORTRAN is defined (which is usually the case),
/// the raw data will be stored with base 1 for Fortran interoperability.
/// The operator(row,col) will then do the conversion to base 0.
///
/// In the first mode of construction, the connectivity table cannot be resized.
/// In the second mode of construction, resizing is possible
class BlockConnectivityImpl {
private:
    friend class IrregularConnectivityImpl;
    friend class MultiBlockConnectivityImpl;
    BlockConnectivityImpl(idx_t rows, idx_t cols, idx_t values[], bool dummy);

public:
    //-- Constructors

    /// @brief Construct connectivity table that needs resizing a-posteriori
    /// Data is owned
    BlockConnectivityImpl();
    BlockConnectivityImpl(idx_t rows, idx_t cols, const std::initializer_list<idx_t>&);

    /// @brief Construct connectivity table wrapping existing raw data.
    /// No resizing can be performed as data is not owned.
    BlockConnectivityImpl(idx_t rows, idx_t cols, idx_t values[]);

#if ATLAS_HIC_COMPILER
    /// @brief Copy ctr (only to be used when calling a cuda kernel)
    // This ctr has to be defined in the header, since __CUDACC__ will identify
    // whether it is compiled it for a GPU kernel
    BlockConnectivityImpl(const BlockConnectivityImpl& other):
        owns_(false),
        values_(other.values_),
        rows_(other.rows_),
        cols_(other.cols_),
        missing_value_(other.missing_value_) {}
#endif

    BlockConnectivityImpl(BlockConnectivityImpl&& other);

    BlockConnectivityImpl& operator=(const BlockConnectivityImpl& other);
    BlockConnectivityImpl& operator=(BlockConnectivityImpl&& other);

    /// @brief Construct a mesh from a Stream (serialization)
    explicit BlockConnectivityImpl(eckit::Stream&);

    /// @brief Destructor
    ~BlockConnectivityImpl();

    ATLAS_HOST_DEVICE
    idx_t index(idx_t i, idx_t j) const;

    //-- Accessors

    /// @brief Access to connectivity table elements for given row and column
    /// The returned index has base 0 regardless if ATLAS_HAVE_FORTRAN is defined.
    ATLAS_HOST_DEVICE
    idx_t operator()(idx_t row_idx, idx_t col_idx) const;

    /// @brief Number of rows
    ATLAS_HOST_DEVICE
    idx_t rows() const { return rows_; }

    /// @brief Number of columns
    ATLAS_HOST_DEVICE
    idx_t cols() const { return cols_; }

    /// @brief Access to raw data.
    /// Note that the connectivity base is 1 in case ATLAS_HAVE_FORTRAN is
    /// defined.
    ATLAS_HOST_DEVICE
    const idx_t* data() const { return values_.data(); }
    ATLAS_HOST_DEVICE
    idx_t* data() { return values_.data(); }

    ATLAS_HOST_DEVICE
    idx_t missing_value() const { return missing_value_; }

    size_t footprint() const;

    //-- Modifiers

    /// @brief Modify row with given values. Values must be given with base 0
    ATLAS_HOST_DEVICE
    void set(idx_t row_idx, const idx_t column_values[]);

    /// @brief Modify (row,col) with given value. Value must be given with base 0
    ATLAS_HOST_DEVICE
    void set(idx_t row_idx, idx_t col_idx, const idx_t value);

    /// @brief Resize connectivity, and add given rows
    /// @note Can only be used when data is owned.
    void add(idx_t rows, idx_t cols, const idx_t values[], bool fortran_array = false);

    bool owns() const { return owns_; }


    friend std::ostream& operator<<(std::ostream& out, const BlockConnectivityImpl& x) {
        x.print(out);
        return out;
    }


protected:
    /// @brief Wrap existing and set owns_ = false
    void rebuild(idx_t rows, idx_t cols, idx_t values[]);


    /// @brief Serialization to Stream
    void encode(eckit::Stream&) const;

    /// @brief Serialization from Stream
    void decode(eckit::Stream&);

    friend eckit::Stream& operator<<(eckit::Stream& s, const BlockConnectivityImpl& x) {
        x.encode(s);
        return s;
    }

    friend eckit::Stream& operator>>(eckit::Stream& s, BlockConnectivityImpl& x) {
        x.decode(s);
        return s;
    }

    void print(std::ostream&) const;

private:
    bool owns_;
    array::SVector<idx_t> values_;

    idx_t rows_;
    idx_t cols_;
    idx_t missing_value_;
};

using IrregularConnectivity  = ConnectivityInterface<IrregularConnectivityImpl>;
using MultiBlockConnectivity = ConnectivityInterface<MultiBlockConnectivityImpl>;
using BlockConnectivity      = BlockConnectivityImpl;
using Connectivity           = IrregularConnectivity;

// -----------------------------------------------------------------------------------------------------

ATLAS_HOST_DEVICE
inline idx_t IrregularConnectivityImpl::operator()(idx_t row_idx, idx_t col_idx) const {
    assert(counts_[row_idx] > (col_idx));
    return values_[displs_[row_idx] + col_idx] FROM_FORTRAN;
}

inline void IrregularConnectivityImpl::set(idx_t row_idx, const idx_t column_values[]) {
    const idx_t N = counts_[row_idx];
    for (idx_t n = 0; n < N; ++n) {
        values_[displs_[row_idx] + n] = column_values[n] TO_FORTRAN;
    }
}

inline void IrregularConnectivityImpl::set(idx_t row_idx, idx_t col_idx, const idx_t value) {
    assert(col_idx < counts_[row_idx]);
    values_[displs_[row_idx] + col_idx] = value TO_FORTRAN;
}

ATLAS_HOST_DEVICE
inline IrregularConnectivityImpl::Row IrregularConnectivityImpl::row(idx_t row_idx) const {
    return IrregularConnectivityImpl::Row(const_cast<idx_t*>(values_.data()) + displs_(row_idx), counts_(row_idx));
}

// -----------------------------------------------------------------------------------------------------

ATLAS_HOST_DEVICE
inline idx_t MultiBlockConnectivityImpl::operator()(idx_t row_idx, idx_t col_idx) const {
    return IrregularConnectivityImpl::operator()(row_idx, col_idx);
}

ATLAS_HOST_DEVICE
inline idx_t MultiBlockConnectivityImpl::operator()(idx_t block_idx, idx_t block_row_idx, idx_t block_col_idx) const {
    return block(block_idx)(block_row_idx, block_col_idx);
}

// -----------------------------------------------------------------------------------------------------

ATLAS_HOST_DEVICE
inline idx_t BlockConnectivityImpl::operator()(idx_t row_idx, idx_t col_idx) const {
    return values_[index(row_idx, col_idx)] FROM_FORTRAN;
}

ATLAS_HOST_DEVICE
inline void BlockConnectivityImpl::set(idx_t row_idx, const idx_t column_values[]) {
    for (idx_t n = 0; n < cols_; ++n) {
        values_[index(row_idx, n)] = column_values[n] TO_FORTRAN;
    }
}

ATLAS_HOST_DEVICE
inline void BlockConnectivityImpl::set(idx_t row_idx, idx_t col_idx, const idx_t value) {
    values_[index(row_idx, col_idx)] = value TO_FORTRAN;
}

ATLAS_HOST_DEVICE
inline idx_t BlockConnectivityImpl::index(idx_t i, idx_t j) const {
    return i * cols_ + j;
}

// ------------------------------------------------------------------------------------------------------

extern "C" {
Connectivity* atlas__Connectivity__create();
MultiBlockConnectivity* atlas__MultiBlockConnectivity__create();
const char* atlas__Connectivity__name(Connectivity* This);
void atlas__Connectivity__rename(Connectivity* This, const char* name);
void atlas__Connectivity__delete(Connectivity* This);
void atlas__Connectivity__displs(Connectivity* This, idx_t*& displs, idx_t& size);
void atlas__Connectivity__counts(Connectivity* This, idx_t*& counts, idx_t& size);
void atlas__Connectivity__values(Connectivity* This, idx_t*& values, idx_t& size);
idx_t atlas__Connectivity__rows(const Connectivity* This);
void atlas__Connectivity__add_values(Connectivity* This, idx_t rows, idx_t cols, idx_t values[]);
void atlas__Connectivity__add_missing(Connectivity* This, idx_t rows, idx_t cols);
idx_t atlas__Connectivity__missing_value(const Connectivity* This);

idx_t atlas__MultiBlockConnectivity__blocks(const MultiBlockConnectivity* This);
BlockConnectivityImpl* atlas__MultiBlockConnectivity__block(MultiBlockConnectivity* This, idx_t block_idx);

idx_t atlas__BlockConnectivity__rows(const BlockConnectivityImpl* This);
idx_t atlas__BlockConnectivity__cols(const BlockConnectivityImpl* This);
idx_t atlas__BlockConnectivity__missing_value(const BlockConnectivityImpl* This);
void atlas__BlockConnectivity__data(BlockConnectivityImpl* This, idx_t*& data, idx_t& rows, idx_t& cols);
void atlas__BlockConnectivity__delete(BlockConnectivityImpl* This);
}

#undef FROM_FORTRAN
#undef TO_FORTRAN
#undef INDEX_DEREF

//------------------------------------------------------------------------------------------------------

}  // namespace mesh

namespace array {

// TODO HACK
// template<> inline DataType::kind_t
// DataType::kind<mesh::BlockConnectivityImpl*>()   { return KIND_INT32;    }
}

}  // namespace atlas
