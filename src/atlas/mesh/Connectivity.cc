/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/Connectivity.h"
#include <limits>
#include "atlas/array.h"
#include "atlas/array/DataType.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/Vector.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/library/defines.h"

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

IrregularConnectivityImpl::IrregularConnectivityImpl( const std::string& name ) :
    owns_( true ),
    data_{array::Array::create<idx_t>( 0 ),    // values
          array::Array::create<size_t>( 1 ),   // displs
          array::Array::create<size_t>( 1 )},  // counts
    values_view_( array::make_host_view<idx_t, 1>( *( data_[_values_] ) ) ),
    displs_view_( array::make_host_view<size_t, 1>( *( data_[_displs_] ) ) ),
    counts_view_( array::make_host_view<size_t, 1>( *( data_[_counts_] ) ) ),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    rows_( 0 ),
    maxcols_( 0 ),
    mincols_( std::numeric_limits<size_t>::max() ),
    ctxt_( 0 ),
    callback_update_( 0 ),
    callback_delete_( 0 ),
    gpu_clone_( this ) {
    rename( name );
    displs_view_( 0 ) = 0;
    counts_view_( 0 ) = 0;
}

// -----------------------------------------------------------------------------

size_t get_total_size_counts( size_t rows, size_t counts[] ) {
    size_t total_size = 0;
    for ( size_t j = 0; j < rows; ++j ) {
        total_size += counts[j];
    }
    return total_size;
}

// -----------------------------------------------------------------------------

IrregularConnectivityImpl::IrregularConnectivityImpl( idx_t values[], size_t rows, size_t displs[], size_t counts[] ) :
    owns_( false ),
    data_{array::Array::wrap<idx_t>( values, array::ArrayShape{get_total_size_counts( rows, counts )} ),
          array::Array::wrap<size_t>( displs, array::ArrayShape{rows} ),
          array::Array::wrap<size_t>( counts, array::ArrayShape{rows} )},
    values_view_( array::make_view<idx_t, 1>( *( data_[_values_] ) ) ),
    displs_view_( array::make_view<size_t, 1>( *( data_[_displs_] ) ) ),
    counts_view_( array::make_view<size_t, 1>( *( data_[_counts_] ) ) ),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    rows_( rows ),
    ctxt_( 0 ),
    callback_update_( 0 ),
    callback_delete_( 0 ),
    gpu_clone_( this ) {
    maxcols_ = 0;
    mincols_ = std::numeric_limits<size_t>::max();
    for ( size_t j = 0; j < rows; ++j ) {
        maxcols_ = std::max( maxcols_, counts[j] );
        mincols_ = std::min( mincols_, counts[j] );
    }
}

IrregularConnectivityImpl::IrregularConnectivityImpl( const IrregularConnectivityImpl& other ) :
    owns_( false ),
    data_{other.data_[0], other.data_[1], other.data_[2]},
    values_view_( other.values_view_ ),
    displs_view_( other.displs_view_ ),
    counts_view_( other.counts_view_ ),
    missing_value_( other.missing_value_ ),
    rows_( other.rows_ ),
    maxcols_( other.maxcols_ ),
    mincols_( other.mincols_ ),
    ctxt_( 0 ),
    gpu_clone_( this ) {}

//------------------------------------------------------------------------------------------------------

IrregularConnectivityImpl::~IrregularConnectivityImpl() {
    on_delete();

    if ( owns_ ) {
        std::for_each( data_.begin(), data_.end(), []( array::Array* a ) {
            assert( a );
            delete a;
            a = 0;
        } );
    }
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::clear() {
    if ( owns() ) {
        data_[_values_]->resize( 0 );
        data_[_displs_]->resize( 1 );
        data_[_counts_]->resize( 1 );
        displs_view_      = array::make_view<size_t, 1>( *( data_[_displs_] ) );
        counts_view_      = array::make_view<size_t, 1>( *( data_[_counts_] ) );
        displs_view_( 0 ) = 0;
        counts_view_( 0 ) = 0;
    }
    else {
        data_[_values_] = nullptr;
        data_[_displs_] = nullptr;
        data_[_counts_] = nullptr;
        // std::for_each(data_.begin(), data_.end(), [](array::Array* a){ a=0;});
    }
    maxcols_ = 0;
    mincols_ = std::numeric_limits<size_t>::max();
    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::on_delete() {
    if ( ctxt_ && callback_delete_ ) callback_delete_( ctxt_ );
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::on_update() {
    if ( ctxt_ && callback_update_ ) callback_update_( ctxt_ );
}

void IrregularConnectivityImpl::resize( size_t old_size, size_t new_size, bool initialize, const idx_t values[],
                                        bool fortran_array ) {
    data_[_values_]->resize( new_size );
    values_view_ = array::make_view<idx_t, 1>( *( data_[_values_] ) );
    // TODO WILLEM isnt this if the other way around?
    idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
    if ( initialize ) {
        for ( size_t j = 0, c = old_size; c < new_size; ++c, ++j ) {
            values_view_( c ) = values[j] + add_base;
        }
    }
    else {
        for ( size_t j = old_size; j < new_size; ++j ) {
            values_view_( j ) = missing_value() TO_FORTRAN;
        }
    }
}

//------------------------------------------------------------------------------------------------------
void IrregularConnectivityImpl::add( size_t rows, size_t cols, const idx_t values[], bool fortran_array ) {
    if ( !owns_ ) throw eckit::AssertionFailed( "HybridConnectivity must be owned to be resized directly" );
    size_t old_size = data_[_values_]->size();

    if ( rows_ == 0 ) old_size = 0;

    size_t new_size = old_size + rows * cols;
    size_t new_rows = rows_ + rows;

    ASSERT( data_[_displs_] != 0 );
    ASSERT( data_[_counts_] != 0 );
    data_[_displs_]->resize( new_rows + 1 );
    data_[_counts_]->resize( new_rows + 1 );
    displs_view_ = array::make_view<size_t, 1>( *( data_[_displs_] ) );
    counts_view_ = array::make_view<size_t, 1>( *( data_[_counts_] ) );

    for ( size_t j = 0; rows_ < new_rows; ++rows_, ++j ) {
        displs_view_( rows_ + 1 ) = displs_view_( rows_ ) + cols;
        counts_view_( rows_ )     = cols;
    }

    maxcols_ = std::max( maxcols_, cols );
    mincols_ = std::min( mincols_, cols );

    resize( old_size, new_size, true, values, fortran_array );

    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::add( const BlockConnectivityImpl& block ) {
    if ( !owns_ ) throw eckit::AssertionFailed( "HybridConnectivity must be owned to be resized directly" );
    bool fortran_array  = FORTRAN_BASE;
    const size_t rows   = block.rows();
    const size_t cols   = block.cols();
    const idx_t* values = block.data();

    std::vector<idx_t> values_vector;
    if ( !block.values_view_.contiguous() ) {
        values_vector.resize( rows * cols );
        values = values_vector.data();
        for ( int i = 0, c = 0; i < rows; ++i ) {
            for ( int j = 0; j < cols; ++j ) {
                values_vector[c++] = block( i, j );
            }
        }
        fortran_array = false;
    }

    add( rows, cols, values, fortran_array );
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::add( size_t rows, const size_t cols[] ) {
    if ( !owns_ ) throw eckit::AssertionFailed( "HybridConnectivity must be owned to be resized directly" );
    size_t old_size = data_[_values_]->size();
    size_t new_size = old_size;
    for ( size_t j = 0; j < rows; ++j )
        new_size += cols[j];
    size_t new_rows = rows_ + rows;
    data_[_displs_]->resize( new_rows + 1 );
    data_[_counts_]->resize( new_rows + 1 );
    displs_view_ = array::make_view<size_t, 1>( *( data_[_displs_] ) );
    counts_view_ = array::make_view<size_t, 1>( *( data_[_counts_] ) );

    for ( size_t j = 0; rows_ < new_rows; ++rows_, ++j ) {
        // TODO isnt this a bug ? I dont understand
        displs_view_( rows_ + 1 ) = displs_view_( rows_ ) + cols[j];
        counts_view_( rows_ )     = cols[j];
        maxcols_                  = std::max( maxcols_, cols[j] );
        mincols_                  = std::min( mincols_, cols[j] );
    }

    resize( old_size, new_size, false, NULL, false );

    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::add( size_t rows, size_t cols ) {
    if ( !owns_ ) throw eckit::AssertionFailed( "HybridConnectivity must be owned to be resized directly" );
    size_t old_size = data_[_values_]->size();

    if ( rows_ == 0 ) old_size = 0;

    size_t new_size = old_size + rows * cols;
    size_t new_rows = rows_ + rows;

    ASSERT( data_[_displs_] != 0 );
    ASSERT( data_[_counts_] != 0 );
    data_[_displs_]->resize( new_rows + 1 );
    data_[_counts_]->resize( new_rows + 1 );
    displs_view_ = array::make_view<size_t, 1>( *( data_[_displs_] ) );
    counts_view_ = array::make_view<size_t, 1>( *( data_[_counts_] ) );

    for ( size_t j = 0; rows_ < new_rows; ++rows_, ++j ) {
        displs_view_( rows_ + 1 ) = displs_view_( rows_ ) + cols;
        counts_view_( rows_ )     = cols;
    }

    maxcols_ = std::max( maxcols_, cols );
    mincols_ = std::min( mincols_, cols );

    const bool dummy_arg_fortran_array = false;
    const idx_t* dummy_arg_values      = NULL;
    resize( old_size, new_size, false, dummy_arg_values, dummy_arg_fortran_array );

    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::insert( size_t position, size_t rows, size_t cols, const idx_t values[],
                                        bool fortran_array ) {
    if ( !owns_ ) throw eckit::AssertionFailed( "HybridConnectivity must be owned to be resized directly" );
    size_t position_displs = displs_view_( position );
    data_[_displs_]->insert( position, rows );
    data_[_counts_]->insert( position, rows );
    displs_view_ = array::make_view<size_t, 1>( *( data_[_displs_] ) );
    counts_view_ = array::make_view<size_t, 1>( *( data_[_counts_] ) );

    displs_view_( position ) = position_displs;
    for ( size_t jrow = position; jrow < position + rows; ++jrow ) {
        counts_view_( jrow ) = cols;
    }
    for ( size_t jrow = position; jrow < displs_view_.size() - 1; ++jrow ) {
        displs_view_( jrow + 1 ) = displs_view_( jrow ) + counts_view_( jrow );
    }
    maxcols_ = std::max( maxcols_, cols );
    mincols_ = std::min( mincols_, cols );

    data_[_values_]->insert( position_displs, rows * cols );
    values_view_ = array::make_view<idx_t, 1>( *( data_[_values_] ) );

    if ( values == NULL ) {
        for ( size_t c = position_displs; c < position_displs + rows * cols; ++c ) {
            values_view_( c ) = missing_value() TO_FORTRAN;
        }
    }
    else {
        idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
        for ( size_t c = position_displs; c < position_displs + rows * cols; ++c ) {
            values_view_( c ) = values[c - position_displs] + add_base;
        }
    }
    rows_ += rows;

    on_update();
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::insert( size_t position, size_t rows, size_t cols ) {
    IrregularConnectivityImpl::insert( position, rows, cols, NULL, false );
}

//------------------------------------------------------------------------------------------------------

void IrregularConnectivityImpl::insert( size_t position, size_t rows, const size_t cols[] ) {
    if ( !owns_ ) throw eckit::AssertionFailed( "HybridConnectivity must be owned to be resized directly" );
    size_t position_displs = displs_view_( position );

    if ( rows_ == 0 ) {
        if ( position > 1 ) {
            data_[_displs_]->insert( position - 1, rows );
            data_[_counts_]->insert( position - 1, rows );
        }
    }
    else {
        data_[_displs_]->insert( position, rows );
        data_[_counts_]->insert( position, rows );
    }
    displs_view_ = array::make_view<size_t, 1>( *( data_[_displs_] ) );
    counts_view_ = array::make_view<size_t, 1>( *( data_[_counts_] ) );

    displs_view_( position ) = position_displs;
    for ( size_t jrow = position; jrow < position + rows; ++jrow ) {
        counts_view_( jrow ) = cols[jrow - position];
        maxcols_             = std::max( maxcols_, counts_view_( jrow ) );
        mincols_             = std::min( mincols_, counts_view_( jrow ) );
    }
    for ( size_t jrow = position; jrow < displs_view_.size() - 1; ++jrow ) {
        displs_view_( jrow + 1 ) = displs_view_( jrow ) + counts_view_( jrow );
    }

    size_t insert_size( 0 );
    for ( size_t j = 0; j < rows; ++j )
        insert_size += cols[j];

    data_[_values_]->insert( position_displs, insert_size );
    values_view_ = array::make_view<idx_t, 1>( *( data_[_values_] ) );

    for ( size_t c = position_displs; c < position_displs + insert_size; ++c ) {
        values_view_( c ) = missing_value() TO_FORTRAN;
    }

    rows_ += rows;
    on_update();
}

void IrregularConnectivityImpl::cloneToDevice() {
    std::for_each( data_.begin(), data_.end(), []( array::Array* a ) { a->cloneToDevice(); } );
    values_view_ = array::make_device_view<idx_t, 1>( *( data_[_values_] ) );
    displs_view_ = array::make_device_view<size_t, 1>( *( data_[_displs_] ) );
    counts_view_ = array::make_device_view<size_t, 1>( *( data_[_counts_] ) );
    gpu_clone_.cloneToDevice();
}
void IrregularConnectivityImpl::cloneFromDevice() {
    std::for_each( data_.begin(), data_.end(), []( array::Array* a ) { a->cloneFromDevice(); } );
    values_view_ = array::make_host_view<idx_t, 1>( *( data_[_values_] ) );
    displs_view_ = array::make_host_view<size_t, 1>( *( data_[_displs_] ) );
    counts_view_ = array::make_host_view<size_t, 1>( *( data_[_counts_] ) );
}
void IrregularConnectivityImpl::syncHostDevice() const {
    std::for_each( data_.begin(), data_.end(), []( array::Array* a ) { a->syncHostDevice(); } );
}
bool IrregularConnectivityImpl::valid() const {
    bool res = true;
    std::for_each( data_.begin(), data_.end(), [&]( array::Array* a ) { res &= a->valid(); } );
    return res;
}
bool IrregularConnectivityImpl::hostNeedsUpdate() const {
    bool res = true;
    std::for_each( data_.begin(), data_.end(), [&]( array::Array* a ) { res &= a->hostNeedsUpdate(); } );
    return res;
}
bool IrregularConnectivityImpl::deviceNeedsUpdate() const {
    bool res = true;
    std::for_each( data_.begin(), data_.end(), [&]( array::Array* a ) { res &= a->deviceNeedsUpdate(); } );
    return res;
}

size_t IrregularConnectivityImpl::footprint() const {
    size_t size = sizeof( *this );
    std::for_each( data_.begin(), data_.end(), [&]( array::Array* a ) { size += a->footprint(); } );
    return size;
}

void IrregularConnectivityImpl::dump( std::ostream& os ) const {
    array::make_host_view<idx_t, 1>( *( data_[_values_] ) ).dump( os );
}

//------------------------------------------------------------------------------------------------------
/*
}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivity::MultiBlockConnectivity( idx_t values[], size_t rows,
size_t displs[], size_t counts[], size_t blocks, size_t block_displs[], size_t
block_cols[] )
  : IrregularConnectivity(values,rows,displs,counts),
    blocks_(blocks),
    block_displs_(array::Array::wrap<size_t>(block_displs,
array::ArrayShape{blocks})),
    block_cols_(array::Array::wrap<size_t>(block_cols,
array::ArrayShape{blocks})),
    block_(blocks),
    block_displs_view_(array::make_view<size_t, 1>(*block_displs_)),
    block_cols_view_(array::make_view<size_t, 1>(*block_cols_))
{
  rebuild_block_connectivity();
}
*/
//------------------------------------------------------------------------------------------------------

MultiBlockConnectivityImpl::MultiBlockConnectivityImpl( const std::string& name ) :
    IrregularConnectivityImpl( name ),
    blocks_( 0 ),
    block_displs_( array::Array::create<size_t>( 1 ) ),
    block_cols_( array::Array::create<size_t>( 1 ) ),
    block_displs_view_( array::make_view<size_t, 1>( *block_displs_ ) ),
    block_cols_view_( array::make_view<size_t, 1>( *block_cols_ ) ),
    block_( 0 ),
    block_view_( make_host_vector_view( block_ ) ),
    gpu_clone_( this ) {
    block_displs_view_( 0 ) = 0;
}

//------------------------------------------------------------------------------------------------------

MultiBlockConnectivityImpl::~MultiBlockConnectivityImpl() {
    if ( block_displs_ ) delete block_displs_;
    if ( block_cols_ ) delete block_cols_;
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::clear() {
    // TODO
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

void MultiBlockConnectivityImpl::cloneToDevice() {
    IrregularConnectivityImpl::cloneToDevice();
    block_displs_->cloneToDevice();
    block_cols_->cloneToDevice();
    block_displs_view_ = array::make_device_view<size_t, 1>( *( block_displs_ ) );
    block_cols_view_   = array::make_device_view<size_t, 1>( *( block_cols_ ) );

    block_.cloneToDevice();
    block_view_ = make_device_vector_view( block_ );

    gpu_clone_.cloneToDevice();
}

void MultiBlockConnectivityImpl::cloneFromDevice() {
    IrregularConnectivityImpl::cloneFromDevice();
    block_displs_->cloneFromDevice();
    block_cols_->cloneFromDevice();
    block_displs_view_ = array::make_host_view<size_t, 1>( *( block_displs_ ) );
    block_cols_view_   = array::make_host_view<size_t, 1>( *( block_cols_ ) );

    block_.cloneFromDevice();
    block_view_ = make_host_vector_view( block_ );
}

void MultiBlockConnectivityImpl::syncHostDevice() const {
    IrregularConnectivityImpl::syncHostDevice();
    block_displs_->syncHostDevice();
    block_cols_->syncHostDevice();
}

bool MultiBlockConnectivityImpl::valid() const {
    return IrregularConnectivityImpl::valid() && block_displs_->valid() && block_cols_->valid();
}

bool MultiBlockConnectivityImpl::hostNeedsUpdate() const {
    return IrregularConnectivityImpl::hostNeedsUpdate() && block_displs_->hostNeedsUpdate() &&
           block_cols_->hostNeedsUpdate();
}

bool MultiBlockConnectivityImpl::deviceNeedsUpdate() const {
    return IrregularConnectivityImpl::deviceNeedsUpdate() && block_displs_->deviceNeedsUpdate() &&
           block_cols_->deviceNeedsUpdate();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::add( size_t rows, size_t cols, const idx_t values[], bool fortran_array ) {
    if ( !owns() ) throw eckit::AssertionFailed( "MultiBlockConnectivity must be owned to be resized directly" );
    size_t old_rows = this->rows();
    IrregularConnectivityImpl::add( rows, cols, values, fortran_array );

    block_displs_->insert( block_displs_->size(), 1 );
    block_cols_->insert( block_cols_->size(), 1 );
    block_displs_view_ = array::make_view<size_t, 1>( *block_displs_ );
    block_cols_view_   = array::make_view<size_t, 1>( *block_cols_ );
    blocks_++;
    block_displs_view_( block_displs_view_.size() - 1 ) = old_rows + rows;
    block_cols_view_( block_cols_view_.size() - 2 )     = cols;

    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::add( const BlockConnectivityImpl& block ) {
    if ( !owns() ) throw eckit::AssertionFailed( "MultiBlockConnectivity must be owned to be resized directly" );
    IrregularConnectivityImpl::add( block );

    block_view_ = make_host_vector_view( block_ );
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::add( size_t rows, size_t cols ) {
    if ( !owns() ) throw eckit::AssertionFailed( "MultiBlockConnectivity must be owned to be resized directly" );
    size_t old_rows = this->rows();
    IrregularConnectivityImpl::add( rows, cols );

    block_displs_->insert( block_displs_->size(), 1 );
    block_cols_->insert( block_cols_->size(), 1 );
    block_displs_view_ = array::make_view<size_t, 1>( *block_displs_ );
    block_cols_view_   = array::make_view<size_t, 1>( *block_cols_ );
    blocks_++;
    block_displs_view_( block_displs_view_.size() - 1 ) = old_rows + rows;
    block_cols_view_( block_cols_view_.size() - 2 )     = cols;

    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::add( size_t rows, const size_t cols[] ) {
    if ( !owns() ) throw eckit::AssertionFailed( "MultiBlockConnectivity must be owned to be resized directly" );
    size_t min      = std::numeric_limits<size_t>::max();
    size_t max      = 0;
    size_t old_rows = this->rows();

    for ( size_t j = 0; j < rows; ++j ) {
        min = std::min( min, cols[j] );
        max = std::min( max, cols[j] );
    }
    if ( min != max )
        throw eckit::AssertionFailed(
            "MultiBlockConnectivity::add(rows,cols[]): "
            "all elements of cols[] must be identical" );
    IrregularConnectivityImpl::add( rows, cols );

    block_displs_->insert( block_displs_->size(), 1 );
    block_cols_->insert( block_cols_->size(), 1 );
    block_displs_view_ = array::make_view<size_t, 1>( *block_displs_ );
    block_cols_view_   = array::make_view<size_t, 1>( *block_cols_ );
    blocks_++;
    block_displs_view_( block_displs_view_.size() - 1 ) = old_rows;
    block_cols_view_( block_cols_view_.size() - 2 )     = max;

    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::insert( size_t position, size_t rows, size_t cols, const idx_t values[],
                                         bool fortran_array ) {
    if ( !owns() ) throw eckit::AssertionFailed( "MultiBlockConnectivity must be owned to be resized directly" );

    ASSERT( blocks_ );

    long blk_idx = blocks_;
    do {
        blk_idx--;
    } while ( blk_idx >= 0l && block_displs_view_( blk_idx ) >= position && cols != block_cols_view_( blk_idx ) );
    ASSERT( blk_idx >= 0l );
    ASSERT( cols == block( blk_idx ).cols() );

    for ( size_t jblk = blk_idx; jblk < blocks_; ++jblk )
        block_displs_view_( jblk + 1 ) += rows;

    IrregularConnectivityImpl::insert( position, rows, cols, values, fortran_array );
    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::insert( size_t position, size_t rows, size_t cols ) {
    if ( !owns() ) throw eckit::AssertionFailed( "MultiBlockConnectivity must be owned to be resized directly" );

    long blk_idx = blocks_;
    do {
        blk_idx--;
    } while ( blk_idx >= 0l && block_displs_view_( blk_idx ) >= position && cols != block_cols_view_( blk_idx ) );

    ASSERT( blk_idx >= 0l );

    IrregularConnectivityImpl::insert( position, rows, cols );

    for ( size_t jblk = blk_idx; jblk < blocks_; ++jblk )
        block_displs_view_( jblk + 1 ) += rows;
    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::insert( size_t position, size_t rows, const size_t cols[] ) {
    if ( !owns() ) throw eckit::AssertionFailed( "MultiBlockConnectivity must be owned to be resized directly" );
    size_t min = std::numeric_limits<size_t>::max();
    size_t max = 0;
    for ( size_t j = 0; j < rows; ++j ) {
        min = std::min( min, cols[j] );
        max = std::min( max, cols[j] );
    }
    if ( min != max )
        throw eckit::AssertionFailed(
            "MultiBlockConnectivity::add(rows,cols[]): "
            "all elements of cls[] must be identical" );

    long blk_idx = blocks_;
    do {
        blk_idx--;
    } while ( blk_idx >= 0l && block_displs_view_( blk_idx ) >= position && max != block_cols_view_( blk_idx ) );

    ASSERT( blk_idx >= 0l );

    IrregularConnectivityImpl::insert( position, rows, cols );

    for ( size_t jblk = blk_idx; jblk < blocks_; ++jblk )
        block_displs_view_( jblk + 1 ) += rows;
    rebuild_block_connectivity();
}

//------------------------------------------------------------------------------------------------------

void MultiBlockConnectivityImpl::rebuild_block_connectivity() {
    block_.resize( blocks_, 0 );
    block_view_ = make_host_vector_view( block_ );

    for ( size_t b = 0; b < blocks_; ++b ) {
        if ( block_view_[b] ) {
            block_view_[b]->rebuild( block_displs_view_( b + 1 ) - block_displs_view_( b ),  // rows
                                     block_cols_view_( b ),                                  // cols
                                     data() + displs( block_displs_view_( b ) ) );
        }
        else {
            block_view_[b] = new BlockConnectivityImpl( block_displs_view_( b + 1 ) - block_displs_view_( b ),  // rows
                                                        block_cols_view_( b ),                                  // cols
                                                        data() + displs( block_displs_view_( b ) ),
                                                        /*own = */ false );
        }
    }
}

//------------------------------------------------------------------------------------------------------

size_t MultiBlockConnectivityImpl::footprint() const {
    size_t size = IrregularConnectivityImpl::footprint();
    size += block_displs_->footprint();
    size += block_cols_->footprint();

    for ( size_t j = 0; j < block_.size(); ++j ) {
        size += block_view_[j]->footprint();
    }
    return size;
}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl() :
    owns_( true ),
    values_( array::Array::create<idx_t>( 1, 1 ) ),
    values_view_( array::make_view<idx_t, 2>( *values_ ) ),
    rows_( 0 ),
    cols_( 0 ),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    gpu_clone_( this ) {}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl( size_t rows, size_t cols, const std::initializer_list<idx_t>& values ) :
    owns_( true ),
    values_( array::Array::create<idx_t>( 1, 1 ) ),
    values_view_( array::make_view<idx_t, 2>( *values_ ) ),
    rows_( rows ),
    cols_( cols ),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    gpu_clone_( this ) {
    delete values_;
    values_        = array::Array::create<idx_t>( rows_, cols_ );
    values_view_   = array::make_view<idx_t, 2>( *values_ );
    idx_t add_base = FORTRAN_BASE;
    auto v         = values.begin();
    for ( size_t i = 0; i < rows_; ++i ) {
        for ( size_t j = 0; j < cols_; ++j ) {
            values_view_( i, j ) = *( v++ ) + add_base;
        }
    }
    ASSERT( v == values.end() );
}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl( size_t rows, size_t cols, idx_t values[] ) :
    owns_( true ),
    values_( array::Array::create<idx_t>( 1, 1 ) ),
    values_view_( array::make_view<idx_t, 2>( *values_ ) ),
    rows_( rows ),
    cols_( cols ),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    gpu_clone_( this ) {
    delete values_;
    values_      = array::Array::create<idx_t>( rows_, cols_ );
    values_view_ = array::make_view<idx_t, 2>( *values_ );
    if ( values_->size() ) {
        idx_t add_base = FORTRAN_BASE;
        idx_t* v       = &values[0];
        for ( size_t i = 0; i < rows_; ++i ) {
            for ( size_t j = 0; j < cols_; ++j ) {
                values_view_( i, j ) = *( v++ ) + add_base;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::BlockConnectivityImpl( size_t rows, size_t cols, idx_t values[], bool dummy ) :
    owns_( false ),
    values_( array::Array::wrap<idx_t>( values, array::ArrayShape{rows, cols} ) ),
    values_view_( array::make_view<idx_t, 2>( *values_ ) ),
    rows_( rows ),
    cols_( cols ),
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    gpu_clone_( this ) {}

//------------------------------------------------------------------------------------------------------

BlockConnectivityImpl::~BlockConnectivityImpl() {
    if ( owns_ ) {
        assert( values_ );
        delete values_;
    }
}

//------------------------------------------------------------------------------------------------------

void BlockConnectivityImpl::rebuild( size_t rows, size_t cols, idx_t values[] ) {
    ASSERT( not owns_ );
    rows_ = rows;
    cols_ = cols;
    assert( values_ );
    delete values_;
    values_      = array::Array::wrap<idx_t>( values, array::ArrayShape{rows, cols} );
    values_view_ = array::make_view<idx_t, 2>( *values_ );
}

//------------------------------------------------------------------------------------------------------

void BlockConnectivityImpl::add( size_t rows, size_t cols, const idx_t values[], bool fortran_array ) {
    if ( !owns_ ) throw eckit::AssertionFailed( "BlockConnectivity must be owned to be resized directly" );
    if ( cols_ != 0 && cols_ != cols )
        throw eckit::AssertionFailed(
            "Cannot add values with different cols than "
            "already existing in BlockConnectivity" );

    values_->resize( rows_ + rows, cols );
    values_view_ = array::make_view<idx_t, 2>( *values_ );

    idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;

    for ( size_t i = 0, i_old = rows_; i < rows; ++i, ++i_old ) {
        for ( size_t j = 0; j < cols; ++j ) {
            values_view_( i_old, j ) = values[i * cols + j] + add_base;
        }
    }

    rows_ += rows;
    cols_ = cols;
}

//------------------------------------------------------------------------------------------------------

size_t BlockConnectivityImpl::footprint() const {
    size_t size = sizeof( *this );
    if ( owns() ) size += values_->footprint();
    return size;
}

void BlockConnectivityImpl::cloneToDevice() {
    values_->cloneToDevice();
    values_view_ = array::make_device_view<idx_t, 2>( *values_ );
    gpu_clone_.cloneToDevice();
}

void BlockConnectivityImpl::cloneFromDevice() {
    values_->cloneFromDevice();
    values_view_ = array::make_host_view<idx_t, 2>( *values_ );
}

bool BlockConnectivityImpl::valid() const {
    return values_->valid();
}

bool BlockConnectivityImpl::hostNeedsUpdate() const {
    return values_->hostNeedsUpdate();
}

bool BlockConnectivityImpl::deviceNeedsUpdate() const {
    return values_->deviceNeedsUpdate();
}

//------------------------------------------------------------------------------------------------------

class ConnectivityPrivateAccess {
private:
    typedef Connectivity::ctxt_t ctxt_t;
    typedef Connectivity::callback_t callback_t;

public:
    ConnectivityPrivateAccess( Connectivity& connectivity ) : connectivity_( connectivity ) {}
    ctxt_t& ctxt() { return connectivity_.ctxt_; }
    callback_t& callback_update() { return connectivity_.callback_update_; }
    callback_t& callback_delete() { return connectivity_.callback_delete_; }

    // TODO : For now return host-view raw data to Fortran, but this should be
    //        reviewed to also possibly return device-view data
    int* values() { return array::make_view<int, 1>( *connectivity_.data_[Connectivity::_values_] ).data(); }
    size_t* displs() { return array::make_view<size_t, 1>( *connectivity_.data_[Connectivity::_displs_] ).data(); }
    size_t* counts() { return array::make_view<size_t, 1>( *connectivity_.data_[Connectivity::_counts_] ).data(); }

    const char* name() { return connectivity_.name_; }

private:
    Connectivity& connectivity_;
};

//------------------------------------------------------------------------------------------------------

extern "C" {
Connectivity* atlas__Connectivity__create() {
    Connectivity* connectivity = 0;
    ATLAS_ERROR_HANDLING( connectivity = new Connectivity(); );
    return connectivity;
}
void atlas__Connectivity__delete( Connectivity* This ) {
    ATLAS_ERROR_HANDLING( delete This );
}

void atlas__connectivity__register_ctxt( Connectivity* This, Connectivity::ctxt_t ctxt ) {
    ConnectivityPrivateAccess access( *This );
    access.ctxt() = ctxt;
}

int atlas__connectivity__ctxt( Connectivity* This, Connectivity::ctxt_t* ctxt ) {
    ConnectivityPrivateAccess access( *This );
    *ctxt = access.ctxt();
    return bool( access.ctxt() );
}

void atlas__connectivity__register_update( Connectivity* This, Connectivity::callback_t callback ) {
    ConnectivityPrivateAccess access( *This );
    access.callback_update() = callback;
}

void atlas__connectivity__register_delete( Connectivity* This, Connectivity::callback_t callback ) {
    ConnectivityPrivateAccess access( *This );
    access.callback_delete() = callback;
}

void atlas__Connectivity__displs( Connectivity* This, size_t*& displs, size_t& size ) {
    ConnectivityPrivateAccess access( *This );
    displs = access.displs();
    size   = This->rows() + 1;
}

void atlas__Connectivity__counts( Connectivity* This, size_t*& counts, size_t& size ) {
    ConnectivityPrivateAccess access( *This );
    counts = access.counts();
    size   = This->rows();
}

void atlas__Connectivity__values( Connectivity* This, int*& values, size_t& size ) {
    ConnectivityPrivateAccess access( *This );
    values = access.values();
    size   = This->rows() ? access.displs()[This->rows()] + 1 : 0;
}

void atlas__Connectivity__add_values( Connectivity* This, size_t rows, size_t cols, int values[] ) {
    This->add( rows, cols, values, true );
}

void atlas__Connectivity__add_missing( Connectivity* This, size_t rows, size_t cols ) {
    This->add( rows, cols );
}

size_t atlas__Connectivity__rows( const Connectivity* This ) {
    return This->rows();
}

int atlas__Connectivity__missing_value( const Connectivity* This ) {
    return This->missing_value() TO_FORTRAN;
}

MultiBlockConnectivity* atlas__MultiBlockConnectivity__create() {
    MultiBlockConnectivity* connectivity = 0;
    ATLAS_ERROR_HANDLING( connectivity = new MultiBlockConnectivity(); );
    return connectivity;
}

size_t atlas__MultiBlockConnectivity__blocks( const MultiBlockConnectivity* This ) {
    return This->blocks();
}

BlockConnectivityImpl* atlas__MultiBlockConnectivity__block( MultiBlockConnectivity* This, size_t block_idx ) {
    ATLAS_ERROR_HANDLING( ASSERT( This != 0 ) );
    BlockConnectivityImpl* block = &This->block( block_idx );
    ASSERT( block != 0 );
    return block;
}

void atlas__BlockConnectivity__delete( BlockConnectivityImpl* This ) {
    ATLAS_ERROR_HANDLING( delete This );
}

size_t atlas__BlockConnectivity__rows( const BlockConnectivityImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This != 0 ) );
    return This->rows();
}

size_t atlas__BlockConnectivity__cols( const BlockConnectivityImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This != 0 ) );
    return This->cols();
}

int atlas__BlockConnectivity__missing_value( const BlockConnectivityImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This != 0 ) );
    return This->missing_value();
}

void atlas__BlockConnectivity__data( BlockConnectivityImpl* This, int*& data, size_t& rows, size_t& cols ) {
    ATLAS_ERROR_HANDLING( ASSERT( This != 0 ) );
    data = This->data();
    rows = This->rows();
    cols = This->cols();
}

const char* atlas__Connectivity__name( Connectivity* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return ConnectivityPrivateAccess( *This ).name(); );
    return 0;
}

void atlas__Connectivity__rename( Connectivity* This, const char* name ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); This->rename( std::string( name ) ); );
}
}

}  // namespace mesh
}  // namespace atlas
