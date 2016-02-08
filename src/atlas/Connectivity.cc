/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <limits>
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/Array.h"
#include "atlas/Connectivity.h"

#ifdef ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define TO_FORTRAN +1
#else
#define FORTRAN_BASE 0
#define TO_FORTRAN
#endif

namespace atlas {

//------------------------------------------------------------------------------------------------------


IrregularConnectivity::IrregularConnectivity( idx_t values[], size_t rows, size_t displs[], size_t counts[] )
  : owns_(false),
    values_(values),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    rows_(rows),
    displs_(displs),
    counts_(counts)
{
  maxcols_ = 0;
  mincols_ = std::numeric_limits<size_t>::max();
  for( size_t j=0; j<rows; ++j ) {
    maxcols_ = std::max(maxcols_,counts[j]);
    mincols_ = std::min(mincols_,counts[j]);
  }
}

void IrregularConnectivity::clear()
{
  if( owns() )
  {
    owned_values_.clear();
    owned_displs_.resize(1); owned_displs_[0]=0ul;
    owned_counts_.resize(1); owned_counts_[0]=0ul;
  }
  values_ = 0;
  rows_ = 0;
  displs_ = 0;
  counts_ = 0;
  maxcols_ = 0;
  mincols_ = std::numeric_limits<size_t>::max();
}

MultiBlockConnectivity::MultiBlockConnectivity( idx_t values[], size_t rows, size_t displs[], size_t counts[], size_t blocks, size_t block_displs[] )
  : IrregularConnectivity(values,rows,displs,counts),
    blocks_(blocks),
    block_displs_(block_displs)
{
  rebuild_block_connectivity();
}

MultiBlockConnectivity::MultiBlockConnectivity() :
  IrregularConnectivity(),
  blocks_(0),
  block_displs_(0),
  owned_block_displs_(1,0ul)
{}

MultiBlockConnectivity::~MultiBlockConnectivity() {}

void MultiBlockConnectivity::clear()
{
  IrregularConnectivity::clear();
  if( owns() )
  {
    owned_block_displs_.resize(1);
    owned_block_displs_[0]=0ul;
  }
  blocks_ = 0;
  block_displs_ = 0;
  block_.clear();
}

void BlockConnectivity::rebuild( size_t rows, size_t cols, idx_t values[] )
{
  rows_ = rows;
  cols_ = cols;
  values_ = values;
}

BlockConnectivity::BlockConnectivity( size_t rows, size_t cols, idx_t values[] )
  : rows_(rows),
    cols_(cols),
    values_(values),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() )
{
}

IrregularConnectivity::IrregularConnectivity() :
  owns_(true),
  values_(0),
  missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
  rows_(0),
  displs_(0),
  counts_(0),
  owned_displs_(1,0ul),
  owned_counts_(1,0ul),
  maxcols_(0),
  mincols_(std::numeric_limits<size_t>::max())
{}

IrregularConnectivity::~IrregularConnectivity() {}

void IrregularConnectivity::add( size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t old_size = owned_values_.size();
  size_t new_size = old_size + rows*cols;
  size_t new_rows = rows_+rows;
  owned_displs_.resize(new_rows+1);
  owned_counts_.resize(new_rows+1);
  for(size_t j=0; rows_<new_rows; ++rows_, ++j) {
    owned_displs_[rows_+1] = owned_displs_[rows_]+cols;
    owned_counts_[rows_] = cols;
  }
  maxcols_ = std::max(maxcols_,cols);
  mincols_ = std::min(mincols_,cols);

  owned_values_.resize(new_size);
  idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
  for(size_t j=0, c=old_size; c<new_size; ++c, ++j) {
    owned_values_[c] = values[j] + add_base;
  }

  values_ = owned_values_.data();
  displs_ = owned_displs_.data();
  counts_ = owned_counts_.data();
}

void IrregularConnectivity::add( const BlockConnectivity& block )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  bool fortran_array = FORTRAN_BASE;
  add(block.rows(),block.cols(),block.data(),fortran_array);
}


void MultiBlockConnectivity::add(size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(rows,cols,values,fortran_array);
  blocks_++;
  owned_block_displs_.push_back(this->rows());
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}

void IrregularConnectivity::add( size_t rows, const size_t cols[] )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t old_size = owned_values_.size();
  size_t new_size = old_size;
  for( size_t j=0; j<rows; ++j )
    new_size += cols[j];
  size_t new_rows = rows_+rows;
  owned_displs_.resize(new_rows+1);
  owned_counts_.resize(new_rows+1);
  for(size_t j=0; rows_<new_rows; ++rows_, ++j) {
    owned_displs_[rows_+1] = owned_displs_[rows_]+cols[j];
    owned_counts_[rows_] = cols[j];
    maxcols_ = std::max(maxcols_,cols[j]);
    mincols_ = std::min(mincols_,cols[j]);
  }

  owned_values_.resize(new_size);
  for( size_t j=old_size; j<new_size; ++j )
    owned_values_[j] = missing_value() TO_FORTRAN;

  values_ = owned_values_.data();
  displs_ = owned_displs_.data();
  counts_ = owned_counts_.data();
}

void IrregularConnectivity::add( size_t rows, size_t cols )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t old_size = owned_values_.size();
  size_t new_size = old_size + rows*cols;
  size_t new_rows = rows_+rows;
  owned_displs_.resize(new_rows+1);
  owned_counts_.resize(new_rows+1);
  for(size_t j=0; rows_<new_rows; ++rows_, ++j) {
    owned_displs_[rows_+1] = owned_displs_[rows_]+cols;
    owned_counts_[rows_] = cols;
  }
  maxcols_ = std::max(maxcols_,cols);
  mincols_ = std::min(mincols_,cols);

  owned_values_.resize(new_size);
  for( size_t j=old_size; j<new_size; ++j )
    owned_values_[j] = missing_value() TO_FORTRAN;

  values_ = owned_values_.data();
  displs_ = owned_displs_.data();
  counts_ = owned_counts_.data();
}


void MultiBlockConnectivity::add( const BlockConnectivity& block )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(block);
  blocks_++;
  owned_block_displs_.push_back(rows());
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}

void MultiBlockConnectivity::add( size_t rows, size_t cols )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(rows,cols);
  blocks_++;
  owned_block_displs_.push_back(this->rows());
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}


void MultiBlockConnectivity::rebuild_block_connectivity()
{
  block_.resize(blocks_);
  for( size_t b=0; b<blocks_; ++b )
  {
    if( block_[b] ) {
      block_[b]->rebuild(
          block_displs_[b+1]-block_displs_[b], // rows
          counts()[block_displs_[b]],          // cols
          data()+displs()[block_displs_[b]]);
    }
    else {
      block_[b].reset(
         new BlockConnectivity(
          block_displs_[b+1]-block_displs_[b], // rows
          counts()[block_displs_[b]],          // cols
          data()+displs()[block_displs_[b]]) );
    }
  }
}

BlockConnectivity::BlockConnectivity() :
  owns_(true), values_(0), rows_(0), cols_(0),
  missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() )
{
}

void BlockConnectivity::add(size_t rows, size_t cols, const idx_t values[], bool fortran_array)
{
  if( !owns_ )
    throw eckit::AssertionFailed("BlockConnectivity must be owned to be resized directly");
  if( cols_!=0 && cols_!=cols)
    throw eckit::AssertionFailed("Cannot add values with different cols than already existing in BlockConnectivity");

  size_t old_size = rows_*cols_;
  size_t new_size = old_size+rows*cols;
  owned_values_.resize(new_size);
  idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
  for( size_t j=0, c=old_size; c<new_size; ++c, ++j ) {
    owned_values_[c] = values[j] + add_base;
  }

  values_=owned_values_.data();
  rows_+=rows;
  cols_=cols;
}


void IrregularConnectivity::insert( size_t position, size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t position_displs = owned_displs_[position];
  owned_displs_.insert( owned_displs_.begin()+position, rows, position_displs );
  owned_counts_.insert( owned_counts_.begin()+position, rows, cols );
  for( size_t jrow=position; jrow<owned_displs_.size()-1; ++jrow ) {
    owned_displs_[jrow+1] = owned_displs_[jrow] + owned_counts_[jrow];
  }
  maxcols_ = std::max(maxcols_,cols);
  mincols_ = std::min(mincols_,cols);

  owned_values_.insert( owned_values_.begin()+position_displs, values,values+rows*cols );

  if( ! fortran_array )
  {
    for(size_t c=position_displs; c<position_displs+rows*cols; ++c) {
      owned_values_[c] += FORTRAN_BASE;
    }
  }
  rows_ += rows;
  values_ = owned_values_.data();
  displs_ = owned_displs_.data();
  counts_ = owned_counts_.data();
}

void IrregularConnectivity::insert( size_t position, size_t rows, size_t cols )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t position_displs = owned_displs_[position];
  owned_displs_.insert( owned_displs_.begin()+position, rows, position_displs );
  owned_counts_.insert( owned_counts_.begin()+position, rows, cols );
  for( size_t jrow=position; jrow<owned_displs_.size()-1; ++jrow ) {
    owned_displs_[jrow+1] = owned_displs_[jrow] + owned_counts_[jrow];
  }
  maxcols_ = std::max(maxcols_,cols);
  mincols_ = std::min(mincols_,cols);

  owned_values_.insert( owned_values_.begin()+position_displs, rows*cols, missing_value() TO_FORTRAN );
  rows_ += rows;
  values_ = owned_values_.data();
  displs_ = owned_displs_.data();
  counts_ = owned_counts_.data();
}

void IrregularConnectivity::insert( size_t position, size_t rows, const size_t cols[] )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t position_displs = owned_displs_[position];
  owned_displs_.insert( owned_displs_.begin()+position, rows, position_displs );
  owned_counts_.insert( owned_counts_.begin()+position, cols, cols+rows );
  for( size_t jrow=position; jrow<owned_displs_.size()-1; ++jrow ) {
    owned_displs_[jrow+1] = owned_displs_[jrow] + owned_counts_[jrow];
    maxcols_ = std::max(maxcols_,owned_counts_[jrow]);
    mincols_ = std::min(mincols_,owned_counts_[jrow]);
  }
  size_t insert_size(0);
  for( size_t j=0; j<rows; ++j )
    insert_size += cols[j];
  owned_values_.insert( owned_values_.begin()+position_displs, insert_size, missing_value() TO_FORTRAN );
  rows_ += rows;
  values_ = owned_values_.data();
  displs_ = owned_displs_.data();
  counts_ = owned_counts_.data();
}


void MultiBlockConnectivity::insert( size_t position, size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  ASSERT( counts()[std::max(position-1ul,0ul)] == cols );

  size_t blk_idx = blocks_;
  do{ blk_idx--; } while( owned_block_displs_[blk_idx] >= position );
  IrregularConnectivity::insert(position,rows,cols,values,fortran_array);

  for( size_t jblk=blk_idx; jblk<blocks_; ++jblk)
    owned_block_displs_[jblk+1] += rows;
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}

void MultiBlockConnectivity::insert( size_t position, size_t rows, size_t cols )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  ASSERT( counts()[std::max(position-1ul,0ul)] == cols );

  size_t blk_idx = blocks_;
  do{ blk_idx--; } while( owned_block_displs_[blk_idx] >= position );

  IrregularConnectivity::insert(position,rows,cols);

  for( size_t jblk=blk_idx; jblk<blocks_; ++jblk)
    owned_block_displs_[jblk+1] += rows;
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}


extern "C" {

}

}  // namespace atlas

