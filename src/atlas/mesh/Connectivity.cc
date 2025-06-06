/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <limits>
#include <ostream>

#include "eckit/io/Buffer.h"
#include "eckit/serialisation/Stream.h"

#include "atlas/array.h"
#include "atlas/array/DataType.h"
#include "atlas/array/MakeView.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/runtime/Exception.h"

#if ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define TO_FORTRAN +1
#else
#define FORTRAN_BASE 0
#define TO_FORTRAN
#endif

namespace atlas {
namespace mesh {
// -----------------------------------------------------------------------------

IrregularConnectivityImpl::IrregularConnectivityImpl(const std::string& name):
    owns_(true),
    values_(0),
    displs_(1),
    counts_(1),
    missing_value_(std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max()),
    rows_(0),
    maxcols_(0),
    mincols_(std::numeric_limits<idx_t>::max()),
    ctxt_(nullptr),
    callback_update_(nullptr),
    callback_delete_(nullptr) {
    rename(name);
    displs_[0] = 0;
    counts_[0] = 0;
}

// -----------------------------------------------------------------------------

idx_t get_total_size_counts(idx_t rows, idx_t counts[]) {
    idx_t total_size = 0;
    for (idx_t j = 0; j < rows; ++j) {
        total_size += counts[j];
    }
    return total_size;
}

// -----------------------------------------------------------------------------

IrregularConnectivityImpl::IrregularConnectivityImpl(idx_t values[], idx_t rows, idx_t displs[], idx_t counts[]):
    owns_(false),
    //TODO need to create clone if pointers are not cuda managed
    values_(values, get_total_size_counts(rows, counts)),
    displs_(displs, rows),
    counts_(counts, rows),
    missing_value_(std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max()),
    rows_(rows),
    ctxt_(nullptr),
    callback_update_(nullptr),
    callback_delete_(nullptr) {
    maxcols_ = 0;
    mincols_ = std::numeric_limits<idx_t>::max();
    for (idx_t j = 0; j < rows; ++j) {
        maxcols_ = std::max(maxcols_, counts[j]);
        mincols_ = std::min(mincols_, counts[j]);
    }
}

//------------------------------------------------------------------------------------------------------

IrregularConnectivityImpl::IrregularConnectivityImpl(eckit::Stream& s) {
    decode_(s);
}

//------------------------------------------------------------------------------------------------------

IrregularConnectivityImpl::~IrregularConnectivityImpl() {
    on_delete();
    //TODO owns is unsed ?
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::clear() {
    //TODO clean this
    if (owns()) {
        values_.resize(0);
        displs_.resize(1);
        counts_.resize(1);
        displs_(0) = 0;
        counts_(0) = 0;
    }
    else {
        //TODO what to do here
        //        data_[_values_] = nullptr;
        //        data_[_displs_] = nullptr;
        //        data_[_counts_] = nullptr;
        // std::for_each(data_.begin(), data_.end(), [](array::Array* a){ a=0;});
    }
    rows_    = 0;
    maxcols_ = 0;
    mincols_ = std::numeric_limits<idx_t>::max();
    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::on_delete() {
    if (ctxt_ && callback_delete_) {
        callback_delete_(ctxt_);
    }
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::on_update() {
    if (ctxt_ && callback_update_) {
        callback_update_(ctxt_);
    }
}

void IrregularConnectivityImpl::resize(idx_t old_size, idx_t new_size, bool initialize, const idx_t values[],
                                       bool fortran_array) {
    values_.resize(new_size);

    idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
    if (initialize) {
        for (idx_t j = 0, c = old_size; c < new_size; ++c, ++j) {
            values_[c] = values[j] + add_base;
        }
    }
    else {
        for (idx_t j = old_size; j < new_size; ++j) {
            values_[j] = missing_value() TO_FORTRAN;
        }
    }
}

//------------------------------------------------------------------------------------------------------
void IrregularConnectivityImpl::add(idx_t rows, idx_t cols, const idx_t values[], bool fortran_array) {
    ATLAS_ASSERT(owns_, "Connectivity must be owned to be resized directly");
    idx_t old_size = values_.size();

    if (rows_ == 0) {
        old_size = 0;
    }

    idx_t new_size = old_size + rows * cols;
    idx_t new_rows = rows_ + rows;

    //TODO what to do here
    //    ATLAS_ASSERT( displs_] != nullptr );
    //    ATLAS_ASSERT( data_[_counts_] != nullptr );
    displs_.resize(new_rows + 1);
    counts_.resize(new_rows);

    for (idx_t j = 0; rows_ < new_rows; ++rows_, ++j) {
        displs_[rows_ + 1] = displs_[rows_] + cols;
        counts_[rows_]     = cols;
    }

    maxcols_ = std::max(maxcols_, cols);
    mincols_ = std::min(mincols_, cols);

    resize(old_size, new_size, true, values, fortran_array);

    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::add(const BlockConnectivityImpl& block) {
    ATLAS_ASSERT(owns_, "Connectivity must be owned to be resized directly");
    bool fortran_array  = FORTRAN_BASE;
    const idx_t rows    = block.rows();
    const idx_t cols    = block.cols();
    const idx_t* values = block.data();

    add(rows, cols, values, fortran_array);
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::add(idx_t rows, const idx_t cols[]) {
    ATLAS_ASSERT(owns_, "Connectivity must be owned to be resized directly");
    idx_t old_size = values_.size();
    idx_t new_size = old_size;
    for (idx_t j = 0; j < rows; ++j) {
        new_size += cols[j];
    }
    idx_t new_rows = rows_ + rows;
    displs_.resize(new_rows + 1);
    counts_.resize(new_rows);

    for (idx_t j = 0; rows_ < new_rows; ++rows_, ++j) {
        displs_[rows_ + 1] = displs_[rows_] + cols[j];
        counts_[rows_]     = cols[j];
        maxcols_           = std::max(maxcols_, cols[j]);
        mincols_           = std::min(mincols_, cols[j]);
    }

    resize(old_size, new_size, false, nullptr, false);

    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::add(idx_t rows, idx_t cols) {
    ATLAS_ASSERT(owns_, "Connectivity must be owned to be resized directly");
    idx_t old_size = values_.size();

    if (rows_ == 0) {
        old_size = 0;
    }

    idx_t new_size = old_size + rows * cols;
    idx_t new_rows = rows_ + rows;

    //TODO
    //    ATLAS_ASSERT( data_[_displs_] != nullptr );
    //    ATLAS_ASSERT( data_[_counts_] != nullptr );
    displs_.resize(new_rows + 1);
    counts_.resize(new_rows);

    for (idx_t j = 0; rows_ < new_rows; ++rows_, ++j) {
        displs_[rows_ + 1] = displs_[rows_] + cols;
        counts_[rows_]     = cols;
    }

    maxcols_ = std::max(maxcols_, cols);
    mincols_ = std::min(mincols_, cols);

    const bool dummy_arg_fortran_array = false;
    const idx_t* dummy_arg_values      = nullptr;
    resize(old_size, new_size, false, dummy_arg_values, dummy_arg_fortran_array);

    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::insert(idx_t position, idx_t rows, idx_t cols, const idx_t values[],
                                       bool fortran_array) {
    ATLAS_ASSERT(owns_, "Connectivity must be owned to be resized directly");
    idx_t position_displs = displs_[position];
    displs_.insert(position, rows);
    counts_.insert(position, rows);

    displs_[position] = position_displs;
    for (idx_t jrow = position; jrow < position + rows; ++jrow) {
        counts_[jrow] = cols;
    }
    for (idx_t jrow = position; jrow < displs_.size() - 1; ++jrow) {
        displs_[jrow + 1] = displs_[jrow] + counts_[jrow];
    }
    maxcols_ = std::max(maxcols_, cols);
    mincols_ = std::min(mincols_, cols);

    values_.insert(position_displs, rows * cols);

    if (values == nullptr) {
        for (idx_t c = position_displs; c < position_displs + rows * cols; ++c) {
            values_[c] = missing_value() TO_FORTRAN;
        }
    }
    else {
        idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
        for (idx_t c = position_displs; c < position_displs + rows * cols; ++c) {
            values_[c] = values[c - position_displs] + add_base;
        }
    }
    rows_ += rows;

    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::insert(idx_t position, idx_t rows, idx_t cols) {
    IrregularConnectivityImpl::insert(position, rows, cols, nullptr, false);
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::insert(idx_t position, idx_t rows, const idx_t cols[]) {
    ATLAS_ASSERT(owns_, "Connectivity must be owned to be resized directly");
    idx_t position_displs = displs_[position];

    if (rows_ == 0) {
        if (position > 1) {
            displs_.insert(position - 1, rows);
            counts_.insert(position - 1, rows);
        }
    }
    else {
        displs_.insert(position, rows);
        counts_.insert(position, rows);
    }

    displs_[position] = position_displs;
    for (idx_t jrow = position; jrow < position + rows; ++jrow) {
        counts_[jrow] = cols[jrow - position];
        maxcols_      = std::max(maxcols_, counts_[jrow]);
        mincols_      = std::min(mincols_, counts_[jrow]);
    }
    for (idx_t jrow = position; jrow < displs_.size() - 1; ++jrow) {
        displs_[jrow + 1] = displs_[jrow] + counts_[jrow];
    }

    idx_t insert_size(0);
    for (idx_t j = 0; j < rows; ++j) {
        insert_size += cols[j];
    }

    values_.insert(position_displs, insert_size);

    for (idx_t c = position_displs; c < position_displs + insert_size; ++c) {
        values_[c] = missing_value() TO_FORTRAN;
    }

    rows_ += rows;
    on_update();
}

size_t IrregularConnectivityImpl::footprint() const {
    size_t size = sizeof(*this);
    size += values_.footprint();
    size += displs_.footprint();
    size += counts_.footprint();
    return size;
}

void IrregularConnectivityImpl::dump(std::ostream& os) const {
    //TODO dump
}

eckit::Stream& operator>>(eckit::Stream& s, array::SVector<idx_t>& x) {
    size_t size;
    s >> size;
    eckit::Buffer buffer(size * sizeof(idx_t));
    s >> buffer;
    idx_t* data = static_cast<idx_t*>(buffer.data());
    idx_t N     = static_cast<idx_t>(size);
    x.resize(N);
    for (idx_t i = 0; i < N; ++i) {
        x(i) = data[i];
    }
    return s;
}

eckit::Stream& operator<<(eckit::Stream& s, const array::SVector<idx_t>& x) {
    size_t size = x.size();
    s << size;
    eckit::Buffer buffer(reinterpret_cast<const char*>(x.data()), size * sizeof(idx_t));
    s << buffer;
    return s;
}

void BlockConnectivityImpl::print(std::ostream& out) const {
    out << "BlockConnectivity:{rows:" << rows_ << ",cols:" << cols_ << ",values:";
    if (values_.size() == 0) {
        out << "null";
    }
    else {
        out << "[";
        for (idx_t i = 0; i < values_.size(); ++i) {
            out << values_[i] << (i < values_.size() - 1 ? "," : "]");
        }
    }
    out << "}";
}


void IrregularConnectivityImpl::encode_(eckit::Stream& s) const {
    s << name();
    s << values_;
    s << displs_;
    s << counts_;
    s << missing_value_;
    s << rows_;
    s << maxcols_;
    s << mincols_;
}


void IrregularConnectivityImpl::decode_(eckit::Stream& s) {
    std::string name;
    s >> name;
    s >> values_;
    s >> displs_;
    s >> counts_;
    s >> missing_value_;
    s >> rows_;
    s >> maxcols_;
    s >> mincols_;
    if (not name.empty()) {
        rename(name);
    }
    ctxt_ = nullptr;
    owns_ = true;
}

//------------------------------------------------------------------------------------------------------
/*
}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivity::MultiBlockConnectivity( idx_t values[], idx_t rows,
idx_t displs[], idx_t counts[], idx_t blocks, idx_t block_displs[], idx_t
block_cols[] )
  : IrregularConnectivity(values,rows,displs,counts),
    blocks_(blocks),
    block_displs_(array::Array::wrap<idx_t>(block_displs,
array::ArrayShape{blocks})),
    block_cols_(array::Array::wrap<idx_t>(block_cols,
array::ArrayShape{blocks})),
    block_(blocks),
    block_displs_view_(array::make_view<idx_t, 1>(*block_displs_)),
    block_cols_view_(array::make_view<idx_t, 1>(*block_cols_))
{
  rebuild_block_connectivity();
}
*/
//------------------------------------------------------------------------------------------------------

MultiBlockConnectivityImpl::MultiBlockConnectivityImpl(const std::string& name):
    IrregularConnectivityImpl(name), blocks_(0), block_displs_(1), block_cols_(1), block_(0) {
    block_displs_(0) = 0;
}

MultiBlockConnectivityImpl::MultiBlockConnectivityImpl(eckit::Stream& s): IrregularConnectivityImpl(s) {
    decode_(s);
}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivityImpl::~MultiBlockConnectivityImpl() {
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::clear() {
    IrregularConnectivityImpl::clear();
    if (owns()) {
        block_displs_.resize(1);
        block_cols_.resize(1);
        block_displs_(0) = 0ul;
    }
    blocks_ = 0;
    block_.clear();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::add(idx_t rows, idx_t cols, const idx_t values[], bool fortran_array) {
    ATLAS_ASSERT(owns(), "MultiBlockConnectivity must be owned to be resized directly");
    idx_t old_rows = this->rows();
    IrregularConnectivityImpl::add(rows, cols, values, fortran_array);

    for (idx_t b = 0; b < blocks_; ++b) {
        ATLAS_ASSERT(block_[b].owns() == false);
    }

    block_displs_.insert(block_displs_.size(), 1);
    block_cols_.insert(block_cols_.size(), 1);

    blocks_++;
    block_displs_[block_displs_.size() - 1] = old_rows + rows;
    block_cols_[block_cols_.size() - 2]     = cols;

    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::add(const BlockConnectivityImpl& block) {
    ATLAS_ASSERT(owns(), "MultiBlockConnectivity must be owned to be resized directly");
    IrregularConnectivityImpl::add(block);
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::add(idx_t rows, idx_t cols) {
    ATLAS_ASSERT(owns(), "MultiBlockConnectivity must be owned to be resized directly");
    idx_t old_rows = this->rows();
    IrregularConnectivityImpl::add(rows, cols);

    for (idx_t b = 0; b < blocks_; ++b) {
        ATLAS_ASSERT(block_[b].owns() == false);
    }


    block_displs_.insert(block_displs_.size(), 1);
    block_cols_.insert(block_cols_.size(), 1);
    blocks_++;

    block_displs_[block_displs_.size() - 1] = old_rows + rows;
    block_cols_[block_cols_.size() - 2]     = cols;

    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::add(idx_t rows, const idx_t cols[]) {
    ATLAS_ASSERT(owns(), "MultiBlockConnectivity must be owned to be resized directly");
    idx_t min      = std::numeric_limits<idx_t>::max();
    idx_t max      = 0;
    idx_t old_rows = this->rows();

    for (idx_t j = 0; j < rows; ++j) {
        min = std::min(min, cols[j]);
        max = std::min(max, cols[j]);
    }
    ATLAS_ASSERT(min == max,
                 "MultiBlockConnectivity::add(rows,cols[]): "
                 "all elements of cols[] must be identical");
    IrregularConnectivityImpl::add(rows, cols);

    for (idx_t b = 0; b < blocks_; ++b) {
        ATLAS_ASSERT(block_[b].owns() == false);
    }

    block_displs_.insert(block_displs_.size(), 1);
    block_cols_.insert(block_cols_.size(), 1);
    blocks_++;
    block_displs_(block_displs_.size() - 1) = old_rows;
    block_cols_[block_cols_.size() - 2]     = max;

    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::insert(idx_t position, idx_t rows, idx_t cols, const idx_t values[],
                                        bool fortran_array) {
    ATLAS_ASSERT(owns(), "MultiBlockConnectivity must be owned to be resized directly");

    ATLAS_ASSERT(blocks_);

    long blk_idx = blocks_;
    do {
        blk_idx--;
    } while (blk_idx >= 0l && block_displs_[blk_idx] >= position && cols != block_cols_[blk_idx]);
    ATLAS_ASSERT(blk_idx >= 0l);
    ATLAS_ASSERT(cols == block(blk_idx).cols());

    for (idx_t jblk = blk_idx; jblk < blocks_; ++jblk) {
        block_displs_[jblk + 1] += rows;
    }

    IrregularConnectivityImpl::insert(position, rows, cols, values, fortran_array);
    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::insert(idx_t position, idx_t rows, idx_t cols) {
    ATLAS_ASSERT(owns(), "MultiBlockConnectivity must be owned to be resized directly");

    long blk_idx = blocks_;
    do {
        blk_idx--;
    } while (blk_idx >= 0l && block_displs_[blk_idx] >= position && cols != block_cols_[blk_idx]);

    ATLAS_ASSERT(blk_idx >= 0l);

    for (idx_t b = 0; b < blocks_; ++b) {
        ATLAS_ASSERT(block_[b].owns() == false);
    }

    IrregularConnectivityImpl::insert(position, rows, cols);

    for (idx_t jblk = blk_idx; jblk < blocks_; ++jblk) {
        block_displs_[jblk + 1] += rows;
    }
    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::insert(idx_t position, idx_t rows, const idx_t cols[]) {
    ATLAS_ASSERT(owns(), "MultiBlockConnectivity must be owned to be resized directly");

    idx_t min = std::numeric_limits<idx_t>::max();
    idx_t max = 0;
    for (idx_t j = 0; j < rows; ++j) {
        min = std::min(min, cols[j]);
        max = std::min(max, cols[j]);
    }
    ATLAS_ASSERT(min == max,
                 "MultiBlockConnectivity::add(rows,cols[]): "
                 "all elements of cls[] must be identical");

    long blk_idx = blocks_;
    do {
        blk_idx--;
    } while (blk_idx >= 0l && block_displs_[blk_idx] >= position && max != block_cols_[blk_idx]);

    ATLAS_ASSERT(blk_idx >= 0l);

    IrregularConnectivityImpl::insert(position, rows, cols);

    for (idx_t jblk = blk_idx; jblk < blocks_; ++jblk) {
        block_displs_[jblk + 1] += rows;
    }
    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::rebuild_block_connectivity() {
    block_.resize(blocks_);

    for (idx_t b = 0; b < blocks_; ++b) {
        block_[b].rebuild(block_displs_[b + 1] - block_displs_[b],  // rows
                          block_cols_[b],                           // cols
                          values_.data() + displs(block_displs_[b]));
    }
}

//------------------------------------------------------------------------------------------------------

size_t MultiBlockConnectivityImpl::footprint() const {
    size_t size = IrregularConnectivityImpl::footprint();
    size += block_displs_.footprint();
    size += block_cols_.footprint();

    for (idx_t j = 0; j < block_.size(); ++j) {
        size += block_[j].footprint();
    }
    return size;
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::encode_(eckit::Stream& s) const {
    s << blocks_;
    s << block_displs_;
    s << block_cols_;
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::decode_(eckit::Stream& s) {
    s >> blocks_;
    s >> block_displs_;
    s >> block_cols_;
    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl():
    owns_(true),
    values_(0),
    rows_(0),
    cols_(0),
    missing_value_(std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max()) {}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl(idx_t rows, idx_t cols, const std::initializer_list<idx_t>& values):
    owns_(true),
    values_(rows * cols),
    rows_(rows),
    cols_(cols),
    missing_value_(std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max()) {
    idx_t add_base = FORTRAN_BASE;
    auto v         = values.begin();
    for (idx_t i = 0; i < rows_; ++i) {
        for (idx_t j = 0; j < cols_; ++j) {
            values_[index(i, j)] = *(v++) + add_base;
        }
    }
    ATLAS_ASSERT(v == values.end());
}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl(idx_t rows, idx_t cols, idx_t values[]):
    owns_(true),
    values_(rows * cols),
    rows_(rows),
    cols_(cols),
    missing_value_(std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max()) {
    if (values_.size()) {
        idx_t add_base = FORTRAN_BASE;
        idx_t* v       = &values[0];
        for (idx_t i = 0; i < rows_; ++i) {
            for (idx_t j = 0; j < cols_; ++j) {
                values_[index(i, j)] = *(v++) + add_base;
            }
        }
    }
}

BlockConnectivityImpl::BlockConnectivityImpl(eckit::Stream& s) {
    decode(s);
}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl(idx_t rows, idx_t cols, idx_t values[], bool /*dummy*/):
    owns_(false),
    values_(values, rows * cols),
    rows_(rows),
    cols_(cols),
    missing_value_(std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max()) {}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl(BlockConnectivityImpl&& other):
    owns_(other.owns_),
    values_(std::move(other.values_)),
    rows_(other.rows_),
    cols_(other.cols_),
    missing_value_(other.missing_value_) {
    other.owns_ = false;
    rows_       = 0;
    cols_       = 0;
}

BlockConnectivityImpl& BlockConnectivityImpl::operator=(BlockConnectivityImpl&& other) {
    owns_          = other.owns_;
    values_        = std::move(other.values_);
    rows_          = other.rows_;
    cols_          = other.cols_;
    missing_value_ = other.missing_value_;
    other.owns_    = false;
    other.rows_    = 0;
    other.cols_    = 0;
    return *this;
}


BlockConnectivityImpl::~BlockConnectivityImpl() {
    values_.clear();
}

//------------------------------------------------------------------------------------------------------

void BlockConnectivityImpl::rebuild(idx_t rows, idx_t cols, idx_t values[]) {
    owns_   = false;
    rows_   = rows;
    cols_   = cols;
    values_ = array::SVector<idx_t>(values, rows * cols);
}

//------------------------------------------------------------------------------------------------------

void BlockConnectivityImpl::add(idx_t rows, idx_t cols, const idx_t values[], bool fortran_array) {
    ATLAS_ASSERT(owns(), "BlockConnectivity must be owned to be resized directly");
    if (cols_ != 0 && cols_ != cols) {
        ATLAS_ASSERT(false,
                     "Cannot add values with different cols than "
                     "already existing in BlockConnectivity");
    }

    values_.resize((rows_ + rows) * cols);
    const idx_t oldrows = rows_;
    idx_t add_base      = fortran_array ? 0 : FORTRAN_BASE;

    rows_ += rows;
    cols_ = cols;

    for (idx_t i = 0; i < rows; ++i) {
        for (idx_t j = 0; j < cols; ++j) {
            values_[index(i + oldrows, j)] = values[i * cols + j] + add_base;
        }
    }
}

void BlockConnectivityImpl::encode(eckit::Stream& s) const {
    s << values_;
    s << rows_;
    s << cols_;
    s << missing_value_;
}

void BlockConnectivityImpl::decode(eckit::Stream& s) {
    owns_ = true;
    s >> values_;
    s >> rows_;
    s >> cols_;
    s >> missing_value_;
}

//------------------------------------------------------------------------------------------------------

size_t BlockConnectivityImpl::footprint() const {
    size_t size = sizeof(*this);
    if (owns()) {
        size += values_.footprint();
    }
    return size;
}

//------------------------------------------------------------------------------------------------------

class ConnectivityPrivateAccess {
private:
    using ctxt_t     = Connectivity::ctxt_t;
    using callback_t = Connectivity::callback_t;

public:
    ConnectivityPrivateAccess(Connectivity& connectivity): connectivity_(connectivity) {}
    ctxt_t& ctxt() { return connectivity_.ctxt_; }
    callback_t& callback_update() { return connectivity_.callback_update_; }
    callback_t& callback_delete() { return connectivity_.callback_delete_; }

    idx_t* values() { return connectivity_.values_.data(); }
    idx_t* displs() { return connectivity_.displs_.data(); }
    idx_t* counts() { return connectivity_.counts_.data(); }

    const char* name() { return connectivity_.name_; }

private:
    Connectivity& connectivity_;
};

//------------------------------------------------------------------------------------------------------

extern "C" {
Connectivity* atlas__Connectivity__create() {
    Connectivity* connectivity = nullptr;
    connectivity               = new Connectivity();
    return connectivity;
}
void atlas__Connectivity__delete(Connectivity* This) {
    delete This;
}

void atlas__connectivity__register_ctxt(Connectivity* This, Connectivity::ctxt_t ctxt) {
    ConnectivityPrivateAccess access(*This);
    access.ctxt() = ctxt;
}

int atlas__connectivity__ctxt(Connectivity* This, Connectivity::ctxt_t* ctxt) {
    ConnectivityPrivateAccess access(*This);
    *ctxt = access.ctxt();
    return bool(access.ctxt());
}

void atlas__connectivity__register_update(Connectivity* This, Connectivity::callback_t callback) {
    ConnectivityPrivateAccess access(*This);
    access.callback_update() = callback;
}

void atlas__connectivity__register_delete(Connectivity* This, Connectivity::callback_t callback) {
    ConnectivityPrivateAccess access(*This);
    access.callback_delete() = callback;
}

void atlas__Connectivity__displs(Connectivity* This, idx_t*& displs, idx_t& size) {
    ConnectivityPrivateAccess access(*This);
    displs = access.displs();
    size   = This->rows() + 1;
}

void atlas__Connectivity__counts(Connectivity* This, idx_t*& counts, idx_t& size) {
    ConnectivityPrivateAccess access(*This);
    counts = access.counts();
    size   = This->rows();
}

void atlas__Connectivity__values(Connectivity* This, idx_t*& values, idx_t& size) {
    ConnectivityPrivateAccess access(*This);
    values = access.values();
    size   = This->rows() ? access.displs()[This->rows()] + 1 : 0;
}

void atlas__Connectivity__add_values(Connectivity* This, idx_t rows, idx_t cols, idx_t values[]) {
    This->add(rows, cols, values, true);
}

void atlas__Connectivity__add_missing(Connectivity* This, idx_t rows, idx_t cols) {
    This->add(rows, cols);
}

idx_t atlas__Connectivity__rows(const Connectivity* This) {
    return This->rows();
}

idx_t atlas__Connectivity__missing_value(const Connectivity* This) {
    return This->missing_value() TO_FORTRAN;
}

MultiBlockConnectivity* atlas__MultiBlockConnectivity__create() {
    MultiBlockConnectivity* connectivity = nullptr;
    connectivity                         = new MultiBlockConnectivity();
    return connectivity;
}

idx_t atlas__MultiBlockConnectivity__blocks(const MultiBlockConnectivity* This) {
    return This->blocks();
}

BlockConnectivityImpl* atlas__MultiBlockConnectivity__block(MultiBlockConnectivity* This, idx_t block_idx) {
    ATLAS_ASSERT(This != nullptr);
    BlockConnectivityImpl* block = &This->block(block_idx);
    ATLAS_ASSERT(block != nullptr);
    return block;
}

void atlas__BlockConnectivity__delete(BlockConnectivityImpl* This) {
    delete This;
}

idx_t atlas__BlockConnectivity__rows(const BlockConnectivityImpl* This) {
    ATLAS_ASSERT(This != nullptr);
    return This->rows();
}

idx_t atlas__BlockConnectivity__cols(const BlockConnectivityImpl* This) {
    ATLAS_ASSERT(This != nullptr);
    return This->cols();
}

idx_t atlas__BlockConnectivity__missing_value(const BlockConnectivityImpl* This) {
    ATLAS_ASSERT(This != nullptr);
    return This->missing_value();
}

void atlas__BlockConnectivity__data(BlockConnectivityImpl* This, idx_t*& data, idx_t& rows, idx_t& cols) {
    ATLAS_ASSERT(This != nullptr);
    data = This->data();
    rows = This->rows();
    cols = This->cols();
}

const char* atlas__Connectivity__name(Connectivity* This) {
    ATLAS_ASSERT(This);
    return ConnectivityPrivateAccess(*This).name();
}

void atlas__Connectivity__rename(Connectivity* This, const char* name) {
    ATLAS_ASSERT(This);
    This->rename(std::string(name));
}
}

}  // namespace mesh
}  // namespace atlas
