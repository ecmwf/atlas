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
#include "atlas/util/runtime/ErrorHandling.h"
#include "atlas/array/Array.h"

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
  values_(0),
  missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
  rows_(0),
  displs_(0),
  counts_(0),
  owned_displs_(1,0ul),
  owned_counts_(1,0ul),
  maxcols_(0),
  mincols_(std::numeric_limits<size_t>::max()),
  ctxt_update_(0),
  ctxt_set_(0),
  ctxt_delete_(0),
  callback_update_(0),
  callback_set_(0),
  callback_delete_(0)
{}

// -----------------------------------------------------------------------------

IrregularConnectivity::IrregularConnectivity( idx_t values[], size_t rows, size_t displs[], size_t counts[] )
  : name_(),
    owns_(false),
    values_(values),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    rows_(rows),
    displs_(displs),
    counts_(counts),
    ctxt_update_(0),
    ctxt_set_(0),
    ctxt_delete_(0),
    callback_update_(0),
    callback_set_(0),
    callback_delete_(0)
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
    owned_values_.clear();
    owned_displs_.resize(1); owned_displs_[0]=0ul;
    owned_counts_.resize(1); owned_counts_[0]=0ul;
  }
  values_ = 0;
  rows_   = 0;
  displs_ = 0;
  counts_ = 0;
  maxcols_ = 0;
  mincols_ = std::numeric_limits<size_t>::max();
  on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::on_delete()
{
  if( ctxt_delete_ && callback_delete_ ) callback_delete_(ctxt_delete_);
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivity::on_update()
{
  if( ctxt_update_ && callback_update_ ) callback_update_(ctxt_update_);
}

//------------------------------------------------------------------------------------------------------

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

  on_update();
}

//------------------------------------------------------------------------------------------------------

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

  on_update();
}

//------------------------------------------------------------------------------------------------------

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

  on_update();
}

//------------------------------------------------------------------------------------------------------

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

  on_update();
}

//------------------------------------------------------------------------------------------------------

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

  on_update();
}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivity::MultiBlockConnectivity( idx_t values[], size_t rows, size_t displs[], size_t counts[], size_t blocks, size_t block_displs[] )
  : IrregularConnectivity(values,rows,displs,counts),
    blocks_(blocks),
    block_displs_(block_displs)
{
  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivity::MultiBlockConnectivity(const std::string& name) :
  IrregularConnectivity(name),
  blocks_(0),
  block_displs_(0),
  owned_block_displs_(1,0ul)
{}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivity::~MultiBlockConnectivity() {}

//------------------------------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::add(size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(rows,cols,values,fortran_array);
  blocks_++;
  owned_block_displs_.push_back(this->rows());
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::add( const BlockConnectivity& block )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(block);
  blocks_++;
  owned_block_displs_.push_back(rows());
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivity::add( size_t rows, size_t cols )
{
  if( !owns() ) throw eckit::AssertionFailed("MultiBlockConnectivity must be owned to be resized directly");
  IrregularConnectivity::add(rows,cols);
  blocks_++;
  owned_block_displs_.push_back(this->rows());
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
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
  owned_block_displs_.push_back(this->rows());
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------------------------------

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
  ASSERT( counts()[std::max(position-1ul,0ul)] == max );

  size_t blk_idx = blocks_;
  do{ blk_idx--; } while( owned_block_displs_[blk_idx] >= position );

  IrregularConnectivity::insert(position,rows,cols);

  for( size_t jblk=blk_idx; jblk<blocks_; ++jblk)
    owned_block_displs_[jblk+1] += rows;
  block_displs_ = owned_block_displs_.data();
  rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------------------------------


BlockConnectivity::BlockConnectivity() :
  owns_(true), values_(0), rows_(0), cols_(0),
  missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() )
{
}

//------------------------------------------------------------------------------------------------------

BlockConnectivity::BlockConnectivity( size_t rows, size_t cols, idx_t values[] )
  : rows_(rows),
    cols_(cols),
    values_(values),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() )
{
}

//------------------------------------------------------------------------------------------------------

void BlockConnectivity::rebuild( size_t rows, size_t cols, idx_t values[] )
{
  rows_ = rows;
  cols_ = cols;
  values_ = values;
}

//------------------------------------------------------------------------------------------------------

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
  idx_t       *values()   { return connectivity_.values_; }
  size_t      *displs()   { return connectivity_.displs_; }
  size_t      *counts()   { return connectivity_.counts_; }

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
  ConnectivityPrivateAccess access(*This);
  displs = access.displs();
  size = This->rows()+1;
}

void atlas__Connectivity__counts(Connectivity* This, size_t* &counts, size_t &size)
{
  ConnectivityPrivateAccess access(*This);
  counts = access.counts();
  size = This->rows();
}

void atlas__Connectivity__values(Connectivity* This, int* &values, size_t &size)
{
  ConnectivityPrivateAccess access(*This);
  values = access.values();
  size = This->rows() ? access.displs()[This->rows()]+1 : 0 ;
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

