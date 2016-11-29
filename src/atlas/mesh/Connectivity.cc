/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <limits>
#include "atlas/mesh/Connectivity.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/DataType.h"

#ifdef ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define TO_FORTRAN +1
#else
#define FORTRAN_BASE 0
#define TO_FORTRAN
#endif

namespace atlas {
namespace mesh {
// -----------------------------------------------------------------------------


IrregularConnectivity::IrregularConnectivity(const std::string& name ) :
  name_(name),
  owns_(true),
  data_{array::Array::create<idx_t>(0), array::Array::create<size_t>(1), array::Array::create<size_t>(1)},
  missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
  rows_(0),
  maxcols_(0),
  mincols_(std::numeric_limits<size_t>::max()),
  ctxt_update_(0),
  ctxt_set_(0),
  ctxt_delete_(0),
  callback_update_(0),
  callback_set_(0),
  callback_delete_(0),
  values_view_(array::make_view<idx_t, 1>(*(data_[_values_]))),
  displs_view_(array::make_view<size_t, 1>(*(data_[_displs_]))),
  counts_view_(array::make_view<size_t, 1>(*(data_[_counts_])))
{
    displs_view_(0) = 0;
    counts_view_(0) = 0;
}

GT_FUNCTION
IrregularConnectivity::IrregularConnectivity(const IrregularConnectivity &other) :
    owns_(other.owns_),
    missing_value_(other.missing_value_),
    rows_(other.rows_),
    maxcols_(other.maxcols_),
    mincols_(other.mincols_),
#ifdef __CUDACC__
    data_{0,0,0},
    values_view_(array::make_device_view<idx_t, 1>(*(other.data_[_values_]))),
    displs_view_(array::make_device_view<size_t, 1>(*(other.data_[_displs_]))),
    counts_view_(array::make_device_view<size_t, 1>(*(other.data_[_counts_])))
#else
    data_{other.data_[0], other.data_[1], other.data_[2]},
    values_view_(array::make_host_view<idx_t, 1>(*(other.data_[_values_]))),
    displs_view_(array::make_host_view<size_t, 1>(*(other.data_[_displs_]))),
    counts_view_(array::make_host_view<size_t, 1>(*(other.data_[_counts_])))
#endif
{}

// -----------------------------------------------------------------------------

size_t get_total_size_counts(size_t rows, size_t counts[])
{
    size_t total_size = 0;
    for( size_t j=0; j<rows; ++j ) {
      total_size += counts[j];
    }
    return total_size;
}

// -----------------------------------------------------------------------------

IrregularConnectivity::IrregularConnectivity( idx_t values[], size_t rows, size_t displs[], size_t counts[] )
  : name_(),
    owns_(false),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    rows_(rows),
    data_{array::Array::wrap<idx_t>(values, array::ArrayShape{get_total_size_counts(rows, counts)}),
            array::Array::wrap<size_t>(displs, array::ArrayShape{rows}),
            array::Array::wrap<size_t>(counts, array::ArrayShape{rows})},
    ctxt_update_(0),
    ctxt_set_(0),
    ctxt_delete_(0),
    callback_update_(0),
    callback_set_(0),
    callback_delete_(0),
    values_view_(array::make_view<idx_t, 1>(*data_[_values_])),
    displs_view_(array::make_view<size_t, 1>(*data_[_displs_])),
    counts_view_(array::make_view<size_t, 1>(*data_[_counts_]))
{
  maxcols_ = 0;
  mincols_ = std::numeric_limits<size_t>::max();
  for( size_t j=0; j<rows; ++j ) {
    maxcols_ = std::max(maxcols_,counts[j]);
    mincols_ = std::min(mincols_,counts[j]);
  }
}

//------------------------------------------------------------------------------------------------------

IrregularConnectivity::~IrregularConnectivity()
{
  on_delete();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::clear()
{
  if( owns() )
  {
      std::for_each(data_.begin(), data_.end(), [](array::Array* a){ assert(a); delete a;});
  }

  std::for_each(data_.begin(), data_.end(), [](array::Array* a){ a=0;});

  maxcols_ = 0;
  mincols_ = std::numeric_limits<size_t>::max();
  on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::on_delete()
{
  if( ctxt_delete_ && callback_delete_ ) callback_delete_(ctxt_delete_);

  if(owns_) {
      std::for_each(data_.begin(), data_.end(), [](array::Array* a){ assert(a); delete a; a = 0;});
  }
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::on_update()
{
  if( ctxt_update_ && callback_update_ ) callback_update_(ctxt_update_);
}

void IrregularConnectivity::resize( size_t old_size, size_t new_size, bool initialize, const idx_t values[], bool fortran_array)
{
  data_[_values_]->resize(new_size);
  values_view_ = array::make_view<idx_t, 1>(*(data_[_values_]));
  //TODO WILLEM isnt this if the other way around?
  idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
  if (initialize) {
    for (size_t j = 0, c = old_size; c < new_size; ++c, ++j) {
      values_view_(c) = values[j] + add_base;
    }
  } else {
    for (size_t j = old_size; j < new_size; ++j) values_view_(j) = missing_value() TO_FORTRAN;
  }
}

//------------------------------------------------------------------------------------------------------
void IrregularConnectivity::add( size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t old_size = data_[_values_]->size();
  size_t new_size = old_size + rows*cols;
  size_t new_rows = rows_+rows;
  data_[_displs_]->resize(new_rows+1);
  data_[_counts_]->resize(new_rows+1);
  displs_view_ = array::make_view<size_t, 1>(*(data_[_displs_]));
  counts_view_ = array::make_view<size_t, 1>(*(data_[_counts_]));

  for(size_t j=0; rows_<new_rows; ++rows_, ++j) {
    displs_view_(rows_+1) = displs_view_(rows_)+cols;
    counts_view_(rows_) = cols;
  }
  maxcols_ = std::max(maxcols_,cols);
  mincols_ = std::min(mincols_,cols);

  resize(old_size, new_size, true, values, fortran_array);

  on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::add( const BlockConnectivity& block )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  bool fortran_array = FORTRAN_BASE;
  add(block.rows(),block.cols(),block.data(),fortran_array);
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::add( size_t rows, const size_t cols[] )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t old_size = data_[_values_]->size();
  size_t new_size = old_size;
  for( size_t j=0; j<rows; ++j )
    new_size += cols[j];
  size_t new_rows = rows_+rows;
  data_[_displs_]->resize(new_rows+1);
  data_[_counts_]->resize(new_rows+1);
  displs_view_ = array::make_view<size_t, 1>(*(data_[_displs_]));
  counts_view_ = array::make_view<size_t, 1>(*(data_[_counts_]));

  for(size_t j=0; rows_<new_rows; ++rows_, ++j) {
    //TODO isnt this a bug ? I dont understand
    displs_view_(rows_+1) = displs_view_(rows_)+cols[j];
    counts_view_(rows_) = cols[j];
    maxcols_ = std::max(maxcols_,cols[j]);
    mincols_ = std::min(mincols_,cols[j]);
  }

  resize(old_size, new_size, false, NULL, false);

  on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::add( size_t rows, size_t cols )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t old_size = data_[_values_]->size();
  size_t new_size = old_size + rows*cols;
  size_t new_rows = rows_+rows;
  data_[_displs_]->resize(new_rows+1);
  data_[_counts_]->resize(new_rows+1);
  displs_view_ = array::make_view<size_t, 1>(*(data_[_displs_]));
  counts_view_ = array::make_view<size_t, 1>(*(data_[_counts_]));

  for(size_t j=0; rows_<new_rows; ++rows_, ++j) {
    displs_view_(rows_+1) = displs_view_(rows_)+cols;
    counts_view_(rows_) = cols;
  }
  maxcols_ = std::max(maxcols_,cols);
  mincols_ = std::min(mincols_,cols);

  resize(old_size, new_size, false, NULL, false);
  on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::insert( size_t position, size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t position_displs = displs_view_(position);
  data_[_displs_]->insert( position, rows);
  data_[_counts_]->insert( position, rows);
  displs_view_ = array::make_view<size_t, 1>(*(data_[_displs_]));
  counts_view_ = array::make_view<size_t, 1>(*(data_[_counts_]));

  displs_view_(position) = position_displs;
  for( size_t jrow=position; jrow<position+rows; ++jrow ) {
      counts_view_(jrow) = cols;
  }
  for( size_t jrow=position; jrow<displs_view_.size()-1; ++jrow ) {
      displs_view_(jrow+1) = displs_view_(jrow) + counts_view_(jrow);
  }
  maxcols_ = std::max(maxcols_,cols);
  mincols_ = std::min(mincols_,cols);

  data_[_values_]->insert(position_displs,rows*cols);
  values_view_ = array::make_view<idx_t, 1>(*(data_[_values_]));

  //TODO WILLEM values was being ignored in the original code
  if(values == NULL) {
      for(size_t c=position_displs; c<position_displs+rows*cols; ++c) {
         values_view_(c) = missing_value() TO_FORTRAN;
      }
  }
  else
  {
    unsigned int base = (fortran_array) ? FORTRAN_BASE : 0;
    for(size_t c=position_displs; c<position_displs+rows*cols; ++c) {
       values_view_(c) = values[c-position_displs] + base;
    }
  }
  rows_ += rows;

  on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::insert( size_t position, size_t rows, size_t cols )
{
    insert(position, rows, cols, NULL, false);
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::insert( size_t position, size_t rows, const size_t cols[] )
{
    if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
    size_t position_displs = displs_view_(position);
    data_[_displs_]->insert( position, rows);
    data_[_counts_]->insert( position, rows);
    displs_view_ = array::make_view<size_t, 1>(*(data_[_displs_]));
    counts_view_ = array::make_view<size_t, 1>(*(data_[_counts_]));

    displs_view_(position) = position_displs;
    for( size_t jrow=position; jrow<position+rows; ++jrow ) {
        counts_view_(jrow) = cols[jrow-position];
        maxcols_ = std::max(maxcols_,counts_view_(jrow));
        mincols_ = std::min(mincols_,counts_view_(jrow));
    }
    for( size_t jrow=position; jrow<displs_view_.size()-1; ++jrow ) {
        displs_view_(jrow+1) = displs_view_(jrow) + counts_view_(jrow);
    }

    size_t insert_size(0);
    for( size_t j=0; j<rows; ++j )
      insert_size += cols[j];

    data_[_values_]->insert(position_displs,insert_size);
    values_view_ = array::make_view<idx_t, 1>(*(data_[_values_]));

    for(size_t c=position_displs; c<position_displs+insert_size; ++c) {
       values_view_(c) = missing_value() TO_FORTRAN;
    }

    rows_ += rows;
    on_update();
}

void IrregularConnectivity::clone_to_device() {
    std::for_each(data_.begin(), data_.end(), [](array::Array* a){ a->clone_to_device();});
    values_view_ = array::make_device_view<idx_t, 1>(*(data_[_values_]));
    displs_view_ = array::make_device_view<size_t, 1>(*(data_[_displs_]));
    counts_view_ = array::make_device_view<size_t, 1>(*(data_[_counts_]));
}
void IrregularConnectivity::clone_from_device() {
    std::for_each(data_.begin(), data_.end(), [](array::Array* a){ a->clone_from_device();});
    values_view_ = array::make_host_view<idx_t, 1>(*(data_[_values_]));
    displs_view_ = array::make_host_view<size_t, 1>(*(data_[_displs_]));
    counts_view_ = array::make_host_view<size_t, 1>(*(data_[_counts_]));
}
bool IrregularConnectivity::valid() const {
    bool res = true;
    std::for_each(data_.begin(), data_.end(), [&](array::Array* a){ res &= a->valid();});
    return res;
}
bool IrregularConnectivity::is_on_host() const {
    bool res=true;
    std::for_each(data_.begin(), data_.end(), [&](array::Array* a){ res &= a->is_on_host();});
    return res;
}
bool IrregularConnectivity::is_on_device() const {
    bool res=true;
    std::for_each(data_.begin(), data_.end(), [&](array::Array* a){ res &= a->is_on_device();});
    return res;
}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivity::MultiBlockConnectivity( idx_t values[], size_t rows, size_t displs[], size_t counts[], size_t blocks, size_t block_displs[], size_t block_cols[] )
  : IrregularConnectivity(values,rows,displs,counts),
    blocks_(blocks),
    block_displs_(array::Array::wrap<size_t>(block_displs, array::ArrayShape{blocks})),
    block_cols_(array::Array::wrap<size_t>(block_cols, array::ArrayShape{blocks})),
    block_(array::Array::create<BlockConnectivity*>(1)),
    block_displs_view_(array::make_view<size_t, 1>(*block_displs_)),
    block_cols_view_(array::make_view<size_t, 1>(*block_cols_)),
    block_view_(array::make_view<BlockConnectivity*, 1>(*block_))
{
//TODO
//  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivity::MultiBlockConnectivity(const std::string& name) :
  IrregularConnectivity(name),
  block_displs_(array::Array::create<size_t>(1)),
  block_cols_(array::Array::create<size_t>(0)),
  block_(array::Array::create<BlockConnectivity*>(0)),
  block_displs_view_(array::make_view<size_t, 1>(*block_displs_)),
  block_cols_view_(array::make_view<size_t, 1>(*block_cols_)),
  block_view_(array::make_view<BlockConnectivity*, 1>(*block_))
{
    block_displs_view_(0) = 0;

}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivity::~MultiBlockConnectivity() {}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::clear()
{
    //TODO
//  IrregularConnectivity::clear();
//  if( owns() )
//  {
//    block_displs_.resize(1);
//    block_displs_[0]=0ul;
//    block_cols_.clear();
//  }
//  blocks_ = 0;
//  block_displs_ = 0;
//  block_cols_ = 0;
//  block_.clear();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::add(size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(rows,cols,values,fortran_array);
  blocks_++;

  block_displs_->insert( block_displs_->size(), 1);
  block_cols_->insert( block_cols_->size(), 1);
  block_displs_view_ = array::make_view<size_t, 1>(*block_displs_);
  block_cols_view_ = array::make_view<size_t, 1>(*block_cols_);

  block_displs_view_(block_displs_view_.size()-1) = this->rows();
  block_cols_view_(block_cols_view_.size()-1) = cols;

//  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::add( const BlockConnectivity& block )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(block);
  blocks_++;

  block_displs_->insert( block_displs_->size(), 1);
  block_cols_->insert( block_cols_->size(), 1);
  block_displs_view_ = array::make_view<size_t, 1>(*block_displs_);
  block_cols_view_ = array::make_view<size_t, 1>(*block_cols_);

  block_displs_view_(block_displs_view_.size()-1) = rows();
  block_cols_view_(block_cols_view_.size()-1) = block.cols();

  //TODO
//  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::add( size_t rows, size_t cols )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(rows,cols);
  blocks_++;

  block_displs_->insert( block_displs_->size(), 1);
  block_cols_->insert( block_cols_->size(), 1);
  block_displs_view_ = array::make_view<size_t, 1>(*block_displs_);
  block_cols_view_ = array::make_view<size_t, 1>(*block_cols_);

  block_displs_view_(block_displs_view_.size()-1) = this->rows();
  block_cols_view_(block_cols_view_.size()-1) = cols;
//TODO
//  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::add( size_t rows, const size_t cols[] )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  size_t min=std::numeric_limits<size_t>::max();
  size_t max=0;
  for( size_t j=0; j<rows; ++j )
  {
    min = std::min(min,cols[j]);
    max = std::min(max,cols[j]);
  }
  if( min != max ) throw eckit::AssertionFailed("MultiBlockConnectivity::add(rows,cols[]): all elements of cls[] must be identical");
  IrregularConnectivity::add(rows,cols);
  blocks_++;

  block_displs_->insert( block_displs_->size(), 1);
  block_cols_->insert( block_cols_->size(), 1);
  block_displs_view_ = array::make_view<size_t, 1>(*block_displs_);
  block_cols_view_ = array::make_view<size_t, 1>(*block_cols_);

  block_displs_view_(block_displs_view_.size()-1) = this->rows();
  block_cols_view_(block_cols_view_.size()-1) = max;
//TODO
//  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::insert( size_t position, size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");

  long blk_idx = blocks_;
  do{ blk_idx--; } while(
     blk_idx >= 0l &&
     block_displs_view_(blk_idx) >= position &&
     cols != block_cols_view_(blk_idx)
  );

  ASSERT( blk_idx >= 0l );

  IrregularConnectivity::insert(position,rows,cols,values,fortran_array);

  for( size_t jblk=blk_idx; jblk<blocks_; ++jblk)
    block_displs_view_(jblk+1) += rows;
// TODO ??
//  block_displs_ = owned_block_displs_.data();
//  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::insert( size_t position, size_t rows, size_t cols )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");

  long blk_idx = blocks_;
  do{ blk_idx--; } while(
     blk_idx >= 0l &&
     block_displs_view_(blk_idx) >= position &&
     cols != block_cols_view_(blk_idx)
  );

  ASSERT( blk_idx >= 0l );

  IrregularConnectivity::insert(position,rows,cols);


  for( size_t jblk=blk_idx; jblk<blocks_; ++jblk)
    block_displs_view_(jblk+1) += rows;
// TODO
//  block_displs_ = owned_block_displs_.data();
//  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::insert( size_t position, size_t rows, const size_t cols[] )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  size_t min=std::numeric_limits<size_t>::max();
  size_t max=0;
  for( size_t j=0; j<rows; ++j )
  {
    min = std::min(min,cols[j]);
    max = std::min(max,cols[j]);
  }
  if( min != max ) throw eckit::AssertionFailed("MultiBlockConnectivity::add(rows,cols[]): all elements of cls[] must be identical");


  long blk_idx = blocks_;
  do{ blk_idx--; } while(
     blk_idx >= 0l &&
     block_displs_view_(blk_idx) >= position &&
     max != block_cols_view_(blk_idx)
  );

  ASSERT( blk_idx >= 0l );

  IrregularConnectivity::insert(position,rows,cols);

  for( size_t jblk=blk_idx; jblk<blocks_; ++jblk)
    block_displs_view_(jblk+1) += rows;
// TODO
//  block_displs_ = owned_block_displs_.data();
//  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::rebuild_block_connectivity()
{
    //TODO
//  block_.resize(blocks_);
//  for( size_t b=0; b<blocks_; ++b )
//  {
//    if( block_[b] ) {
//      block_[b]->rebuild(
//          block_displs_[b+1]-block_displs_[b], // rows
//          block_cols_[b],                      // cols
//          data()+displs()[block_displs_[b]]);
//    }
//    else {
//      block_[b].reset(
//         new BlockConnectivity(
//          block_displs_[b+1]-block_displs_[b], // rows
//          block_cols_[b],                      // cols
//          data()+displs()[block_displs_[b]]) );
//    }
//  }
}

//------------------------------------------------------------------------------------------------------

BlockConnectivity::BlockConnectivity() :
  owns_(true), rows_(0), cols_(0), values_(array::Array::create<idx_t>(0,0) ),
  values_view_(array::make_view<idx_t, 2>(*values_)),
  missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() )
{}

//------------------------------------------------------------------------------------------------------

BlockConnectivity::BlockConnectivity( size_t rows, size_t cols, idx_t values[] )
  : owns_(false),
    rows_(rows),
    cols_(cols),
    values_(array::Array::wrap<idx_t>(values, array::ArrayShape{rows, cols})),
    values_view_(array::make_view<idx_t, 2>(*values_)),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() )
{}

GT_FUNCTION
BlockConnectivity::BlockConnectivity(const BlockConnectivity& other)
    : rows_(other.rows_),
      cols_(other.cols_),
      values_(other.values_),
#ifdef __CUDACC__
      values_view_(array::make_device_view<idx_t, 2>(*values_)),
#else
      values_view_(array::make_device_view<idx_t, 2>(*values_)),
#endif
      missing_value_( other.missing_value_)
{}

//------------------------------------------------------------------------------------------------------

BlockConnectivity::~BlockConnectivity() {
    if(owns_) {
        assert(values_);
        delete values_;
    }
}

//------------------------------------------------------------------------------------------------------

void BlockConnectivity::rebuild( size_t rows, size_t cols, idx_t values[] )
{
  rows_ = rows;
  cols_ = cols;
  assert(values_);
  delete values_;
  values_ = array::Array::wrap<idx_t>(values, array::ArrayShape{rows, cols});
  values_view_ = array::make_view<idx_t, 2>(*values_);
}

//------------------------------------------------------------------------------------------------------

void BlockConnectivity::add(size_t rows, size_t cols, const idx_t values[], bool fortran_array)
{
  if (!owns_) throw eckit::AssertionFailed("BlockConnectivity must be owned to be resized directly");
  if (cols_ != 0 && cols_ != cols)
    throw eckit::AssertionFailed("Cannot add values with different cols than already existing in BlockConnectivity");

  values_->resize(rows_+rows, cols);
  values_view_ = array::make_view<idx_t, 2>(*values_);

  idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;

    for (size_t i = 0, i_old=rows_; i < rows; ++i, ++i_old) {
      for (size_t j = 0; j < cols; ++j) {
        values_view_(i_old, j) = values[i*cols + j] + add_base;
      }
    }

  rows_ += rows;
  cols_ = cols;
}

void BlockConnectivity::clone_to_device()  {
    values_->clone_to_device();
std::cout << "II " << values_view_.data() << std::endl;
    values_view_ = array::make_device_view<idx_t, 2>(*values_);
std::cout << "II " << values_view_.data() << std::endl;

}
void BlockConnectivity::clone_from_device() {
    values_->clone_from_device();
    values_view_ = array::make_host_view<idx_t, 2>(*values_);
}

bool BlockConnectivity::valid() const {
    return values_->valid();
}

bool BlockConnectivity::is_on_host() const {
    return values_->is_on_host();
}

bool BlockConnectivity::is_on_device() const {
    return values_->is_on_device();
}

//------------------------------------------------------------------------------------------------------

class ConnectivityPrivateAccess
{
private:
  typedef Connectivity::ctxt_t     ctxt_t;
  typedef Connectivity::callback_t callback_t;
public:
  ConnectivityPrivateAccess(Connectivity& connectivity) : connectivity_(connectivity)
  {
  }
  ctxt_t      &ctxt_update()     { return connectivity_.ctxt_update_; }
  ctxt_t      &ctxt_set()        { return connectivity_.ctxt_set_; }
  ctxt_t      &ctxt_delete()     { return connectivity_.ctxt_delete_; }
  callback_t  &callback_update() { return connectivity_.callback_update_; }
  callback_t  &callback_set()    { return connectivity_.callback_set_; }
  callback_t  &callback_delete() { return connectivity_.callback_delete_; }
//TODO return array or arrayview?
  array::Array       *values()   { return connectivity_.data_[Connectivity::_values_]; }
  array::Array      *displs()   { return connectivity_.data_[Connectivity::_displs_]; }
  array::Array      *counts()   { return connectivity_.data_[Connectivity::_counts_]; }

private:
  Connectivity& connectivity_;
};

//------------------------------------------------------------------------------------------------------

extern "C"
{
Connectivity* atlas__Connectivity__create()
{
  Connectivity* connectivity = 0;
  ATLAS_ERROR_HANDLING(
    connectivity = new Connectivity();
  );
  return connectivity;
}
void atlas__Connectivity__delete(Connectivity* This)
{
  ATLAS_ERROR_HANDLING(delete This);
}

void atlas__connectivity__register_update(Connectivity* This, Connectivity::callback_t callback, Connectivity::ctxt_t ctxt )
{
  ConnectivityPrivateAccess access(*This);
  access.ctxt_update() = ctxt;
  access.callback_update() = callback;
}

int atlas__connectivity__ctxt_update(Connectivity* This, Connectivity::ctxt_t* ctxt)
{
  ConnectivityPrivateAccess access(*This);
  *ctxt = access.ctxt_update();
  return bool( access.ctxt_update() );
}

void atlas__connectivity__register_delete(Connectivity* This, Connectivity::callback_t callback, Connectivity::ctxt_t ctxt )
{
  ConnectivityPrivateAccess access(*This);
  access.ctxt_delete() = ctxt;
  access.callback_delete() = callback;
}

int atlas__connectivity__ctxt_delete(Connectivity* This, Connectivity::ctxt_t* ctxt)
{
  ConnectivityPrivateAccess access(*This);
  *ctxt = access.ctxt_delete();
  return bool( access.ctxt_delete() );
}

void atlas__Connectivity__displs(Connectivity* This, size_t* &displs, size_t &size)
{
    //TODO
//  ConnectivityPrivateAccess access(*This);
//  displs = access.displs();
//  size = This->rows()+1;
}

void atlas__Connectivity__counts(Connectivity* This, size_t* &counts, size_t &size)
{
    //TODO
//  ConnectivityPrivateAccess access(*This);
//  counts = access.counts();
//  size = This->rows();
}

void atlas__Connectivity__values(Connectivity* This, int* &values, size_t &size)
{
    //TODO
//  ConnectivityPrivateAccess access(*This);
//  values = access.values();
//  size = This->rows() ? access.displs()[This->rows()]+1 : 0 ;
}

void atlas__Connectivity__add_values(Connectivity* This, size_t rows, size_t cols, int values[])
{
  This->add(rows,cols,values,true);
}

void atlas__Connectivity__add_missing(Connectivity* This, size_t rows, size_t cols)
{
  This->add(rows,cols);
}

size_t atlas__Connectivity__rows(const Connectivity* This)
{
  return This->rows();
}

int atlas__Connectivity__missing_value(const Connectivity* This)
{
  return This->missing_value() TO_FORTRAN;
}

MultiBlockConnectivity* atlas__MultiBlockConnectivity__create()
{
  MultiBlockConnectivity* connectivity = 0;
  ATLAS_ERROR_HANDLING(
    connectivity = new MultiBlockConnectivity();
  );
  return connectivity;
}

size_t atlas__MultiBlockConnectivity__blocks(const MultiBlockConnectivity* This)
{
  return This->blocks();
}

BlockConnectivity* atlas__MultiBlockConnectivity__block(MultiBlockConnectivity* This, size_t block_idx)
{
  ATLAS_ERROR_HANDLING( ASSERT(This != 0 ) );
  BlockConnectivity* block = &This->block(block_idx);
  ASSERT( block != 0 );
  return block;
}

void atlas__BlockConnectivity__delete(BlockConnectivity* This)
{
  ATLAS_ERROR_HANDLING( delete This );
}

size_t atlas__BlockConnectivity__rows(const BlockConnectivity* This)
{
  ATLAS_ERROR_HANDLING( ASSERT(This != 0 ) );
  return This->rows();
}

size_t atlas__BlockConnectivity__cols(const BlockConnectivity* This)
{
  ATLAS_ERROR_HANDLING( ASSERT(This != 0 ) );
  return This->cols();
}

int atlas__BlockConnectivity__missing_value(const BlockConnectivity* This)
{
  ATLAS_ERROR_HANDLING( ASSERT(This != 0 ) );
  return This->missing_value();
}

void atlas__BlockConnectivity__data(BlockConnectivity* This, int* &data, size_t &rows, size_t &cols)
{
  ATLAS_ERROR_HANDLING( ASSERT(This != 0 ) );
  data = This->data();
  rows = This->rows();
  cols = This->cols();
}

const char* atlas__Connectivity__name (Connectivity* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->name().c_str();
  );
  return 0;
}

void atlas__Connectivity__rename(Connectivity* This, const char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    This->rename( std::string(name) );
  );
}

}

} // namespace mesh
} // namespace atlas

