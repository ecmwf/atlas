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
}

MultiBlockConnectivity::MultiBlockConnectivity( idx_t values[], size_t rows, size_t displs[], size_t counts[], size_t blocks, size_t block_displs[] )
  : IrregularConnectivity(values,rows,displs,counts),
    blocks_(blocks),
    block_displs_(block_displs)
{
  regenerate_block_connectivity();
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
  owned_counts_(1,0ul)
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
  regenerate_block_connectivity();
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
  regenerate_block_connectivity();
}


void MultiBlockConnectivity::regenerate_block_connectivity()
{
  block_.resize(blocks_);
  for( size_t b=0; b<blocks_; ++b )
  {
    block_[b].reset(
       new BlockConnectivity(
        block_displs_[b+1]-block_displs_[b], // rows
        counts()[block_displs_[b]],          // cols
        data()+displs()[block_displs_[b]]) );
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


extern "C" {

}

}  // namespace atlas

