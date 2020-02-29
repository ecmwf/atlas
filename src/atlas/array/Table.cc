/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Table.h"

#include <algorithm>
#include <limits>

#include "atlas/array.h"
#include "atlas/array/DataType.h"
#include "atlas/library/defines.h"
#include "atlas/runtime/Exception.h"

#if ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define TO_FORTRAN +1
#else
#define FORTRAN_BASE 0
#define TO_FORTRAN
#endif

namespace atlas {
namespace array {

// ----------------------------------------------------------------------------

Table::Table( const std::string& name ) :
    name_( name ),
    owns_( true ),
    data_{Array::create<idx_t>( 0 ),    // values
          Array::create<size_t>( 1 ),   // displs
          Array::create<size_t>( 1 )},  // counts
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    rows_( 0 ),
    maxcols_( 0 ),
    mincols_( std::numeric_limits<size_t>::max() ),
    displs_( make_host_view<size_t, 1>( *( data_[_displs_] ) ) ),
    counts_( make_host_view<size_t, 1>( *( data_[_counts_] ) ) ),
    values_( make_host_view<idx_t, 1>( *( data_[_values_] ) ) ) {
    displs_( 0 ) = 0;
    counts_( 0 ) = 0;
}

// ----------------------------------------------------------------------------

namespace {
static size_t get_total_size_counts( size_t rows, size_t counts[] ) {
    size_t total_size = 0;
    for ( size_t j = 0; j < rows; ++j ) {
        total_size += counts[j];
    }
    return total_size;
}
}  // namespace

// ----------------------------------------------------------------------------

Table::Table( idx_t values[], size_t rows, size_t displs[], size_t counts[] ) :
    name_(),
    owns_( false ),
    data_{Array::wrap<idx_t>( values, array::ArrayShape{get_total_size_counts( rows, counts )} ),
          Array::wrap<size_t>( displs, array::ArrayShape{rows} ),
          Array::wrap<size_t>( counts, array::ArrayShape{rows} )},
    missing_value_( std::numeric_limits<idx_t>::is_signed ? -1 : std::numeric_limits<idx_t>::max() ),
    rows_( rows ),
    maxcols_( 0 ),
    mincols_( std::numeric_limits<size_t>::max() ),
    displs_( make_host_view<size_t, 1>( *( data_[_displs_] ) ) ),
    counts_( make_host_view<size_t, 1>( *( data_[_counts_] ) ) ),
    values_( make_host_view<idx_t, 1>( *( data_[_values_] ) ) ) {
    for ( size_t j = 0; j < rows; ++j ) {
        maxcols_ = std::max( maxcols_, counts[j] );
        mincols_ = std::min( mincols_, counts[j] );
    }
}

// ----------------------------------------------------------------------------

Table::~Table() {
    if ( owns_ ) {
        std::for_each( data_.begin(), data_.end(), []( array::Array* a ) {
            assert( a );
            delete a;
            a = 0;
        } );
    }
}

// ----------------------------------------------------------------------------

void Table::clear() {
    if ( owns() ) {
        resize_values( 0 );
        resize_counts_and_displs( 1 );
        displs_( 0 ) = 0;
        counts_( 0 ) = 0;
    }
    else {
        std::for_each( data_.begin(), data_.end(), []( array::Array* a ) { a = 0; } );
        // data_ and host_views will be invalid now!
    }
    maxcols_ = 0;
    mincols_ = std::numeric_limits<size_t>::max();
}

// ----------------------------------------------------------------------------

void Table::resize_values( size_t old_size, size_t new_size, bool initialize, const idx_t values[],
                           bool fortran_array ) {
    resize_values( new_size );
    idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
    if ( initialize ) {
        for ( size_t j = 0, c = old_size; c < new_size; ++c, ++j ) {
            values_( c ) = values[j] + add_base;
        }
    }
    else {
        for ( size_t j = old_size; j < new_size; ++j ) {
            values_( j ) = missing_value() TO_FORTRAN;
        }
    }
}

// ----------------------------------------------------------------------------

void Table::resize_counts_and_displs( size_t size ) {
    ATLAS_ASSERT( data_[_displs_] != 0 );
    ATLAS_ASSERT( data_[_counts_] != 0 );
    data_[_displs_]->resize( size );
    data_[_counts_]->resize( size );
    displs_ = make_host_view<size_t, 1>( *( data_[_displs_] ) );
    counts_ = make_host_view<size_t, 1>( *( data_[_counts_] ) );
}

// ----------------------------------------------------------------------------

void Table::insert_counts_and_displs( size_t position, size_t rows ) {
    ATLAS_ASSERT( data_[_displs_] != 0 );
    ATLAS_ASSERT( data_[_counts_] != 0 );
    data_[_displs_]->insert( position, rows );
    data_[_counts_]->insert( position, rows );
    displs_ = make_host_view<size_t, 1>( *( data_[_displs_] ) );
    counts_ = make_host_view<size_t, 1>( *( data_[_counts_] ) );
}

// ----------------------------------------------------------------------------

void Table::resize_values( size_t size ) {
    ATLAS_ASSERT( data_[_values_] != 0 );
    data_[_values_]->resize( size );
    values_ = make_host_view<idx_t, 1>( *( data_[_values_] ) );
}

// ----------------------------------------------------------------------------

void Table::insert_values( size_t position, size_t size ) {
    ATLAS_ASSERT( data_[_values_] != 0 );
    data_[_values_]->insert( position, size );
    values_ = make_host_view<idx_t, 1>( *( data_[_values_] ) );
}

// ----------------------------------------------------------------------------

void Table::add( size_t rows, size_t cols, const idx_t values[], bool fortran_array ) {
    ATLAS_ASSERT( owns_, "HybridConnectivity must be owned to be resized directly" );
    size_t old_size = size();

    if ( rows_ == 0 )
        old_size = 0;

    size_t new_size = old_size + rows * cols;
    size_t new_rows = rows_ + rows;

    resize_counts_and_displs( new_rows + 1 );

    for ( size_t j = 0; rows_ < new_rows; ++rows_, ++j ) {
        displs_( rows_ + 1 ) = displs_( rows_ ) + cols;
        counts_( rows_ )     = cols;
    }

    maxcols_ = std::max( maxcols_, cols );
    mincols_ = std::min( mincols_, cols );

    resize_values( old_size, new_size, true, values, fortran_array );
}

// ----------------------------------------------------------------------------

void Table::add( size_t rows, const size_t cols[] ) {
    ATLAS_ASSERT( owns_, "HybridConnectivity must be owned to be resized directly" );
    size_t old_size = size();
    size_t new_size = old_size;
    for ( size_t j = 0; j < rows; ++j )
        new_size += cols[j];
    size_t new_rows = rows_ + rows;
    resize_counts_and_displs( new_rows + 1 );

    for ( size_t j = 0; rows_ < new_rows; ++rows_, ++j ) {
        displs_( rows_ + 1 ) = displs_( rows_ ) + cols[j];
        counts_( rows_ )     = cols[j];
        maxcols_             = std::max( maxcols_, cols[j] );
        mincols_             = std::min( mincols_, cols[j] );
    }

    resize_values( old_size, new_size, false, NULL, false );
}

// ----------------------------------------------------------------------------

void Table::add( size_t rows, size_t cols ) {
    ATLAS_ASSERT( owns_, "HybridConnectivity must be owned to be resized directly" );
    size_t old_size = size();

    if ( rows_ == 0 )
        old_size = 0;

    size_t new_size = old_size + rows * cols;
    size_t new_rows = rows_ + rows;
    resize_counts_and_displs( new_rows + 1 );
    for ( size_t j = 0; rows_ < new_rows; ++rows_, ++j ) {
        displs_( rows_ + 1 ) = displs_( rows_ ) + cols;
        counts_( rows_ )     = cols;
    }

    maxcols_ = std::max( maxcols_, cols );
    mincols_ = std::min( mincols_, cols );

    const bool dummy_arg_fortran_array = false;
    const idx_t* dummy_arg_values      = nullptr;
    resize_values( old_size, new_size, false, dummy_arg_values, dummy_arg_fortran_array );
}

// ----------------------------------------------------------------------------

void Table::insert( size_t position, size_t rows, size_t cols, const idx_t values[], bool fortran_array ) {
    ATLAS_ASSERT( owns_, "HybridConnectivity must be owned to be resized directly" );
    size_t position_displs = displs_( position );
    insert_counts_and_displs( position, rows );

    displs_( position ) = position_displs;
    for ( size_t jrow = position; jrow < position + rows; ++jrow ) {
        counts_( jrow ) = cols;
    }
    for ( size_t jrow = position; jrow < displs_.size() - 1; ++jrow ) {
        displs_( jrow + 1 ) = displs_( jrow ) + counts_( jrow );
    }
    maxcols_ = std::max( maxcols_, cols );
    mincols_ = std::min( mincols_, cols );

    insert_values( position_displs, rows * cols );

    if ( values == nullptr ) {
        for ( size_t c = position_displs; c < position_displs + rows * cols; ++c ) {
            values_( c ) = missing_value() TO_FORTRAN;
        }
    }
    else {
        idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
        for ( size_t c = position_displs; c < position_displs + rows * cols; ++c ) {
            values_( c ) = values[c - position_displs] + add_base;
        }
    }
    rows_ += rows;
}

// ----------------------------------------------------------------------------

void Table::insert( size_t position, size_t rows, size_t cols ) {
    Table::insert( position, rows, cols, nullptr, false );
}

// ----------------------------------------------------------------------------

void Table::insert( size_t position, size_t rows, const size_t cols[] ) {
    ATLAS_ASSERT( owns_, "HybridConnectivity must be owned to be resized directly" );
    size_t position_displs = displs_( position );

    if ( rows_ == 0 ) {
        if ( position > 1 ) {
            insert_counts_and_displs( position - 1, rows );
        }
    }
    else {
        insert_counts_and_displs( position, rows );
    }
    displs_( position ) = position_displs;
    for ( size_t jrow = position; jrow < position + rows; ++jrow ) {
        counts_( jrow ) = cols[jrow - position];
        maxcols_        = std::max( maxcols_, counts_( jrow ) );
        mincols_        = std::min( mincols_, counts_( jrow ) );
    }
    for ( size_t jrow = position; jrow < displs_.size() - 1; ++jrow ) {
        displs_( jrow + 1 ) = displs_( jrow ) + counts_( jrow );
    }

    size_t insert_size( 0 );
    for ( size_t j = 0; j < rows; ++j )
        insert_size += cols[j];

    insert_values( position_displs, insert_size );

    for ( size_t c = position_displs; c < position_displs + insert_size; ++c ) {
        values_( c ) = missing_value() TO_FORTRAN;
    }

    rows_ += rows;
}

// ----------------------------------------------------------------------------

void Table::updateDevice() const {
    std::for_each( data_.begin(), data_.end(), []( array::Array* a ) { a->updateDevice(); } );
}

// ----------------------------------------------------------------------------

void Table::updateHost() const {
    std::for_each( data_.begin(), data_.end(), []( array::Array* a ) { a->updateHost(); } );
}

// ----------------------------------------------------------------------------

void Table::syncHostDevice() const {
    std::for_each( data_.begin(), data_.end(), []( array::Array* a ) { a->syncHostDevice(); } );
}

// ----------------------------------------------------------------------------

bool Table::valid() const {
    bool res = true;
    std::for_each( data_.begin(), data_.end(), [&]( array::Array* a ) { res &= a->valid(); } );
    return res;
}

// ----------------------------------------------------------------------------

bool Table::hostNeedsUpdate() const {
    bool res = true;
    std::for_each( data_.begin(), data_.end(), [&]( array::Array* a ) { res &= a->hostNeedsUpdate(); } );
    return res;
}

// ----------------------------------------------------------------------------

bool Table::deviceNeedsUpdate() const {
    bool res = true;
    std::for_each( data_.begin(), data_.end(), [&]( array::Array* a ) { res &= a->deviceNeedsUpdate(); } );
    return res;
}

// ----------------------------------------------------------------------------

size_t Table::footprint() const {
    size_t size = sizeof( *this );
    if ( owns_ ) {
        std::for_each( data_.begin(), data_.end(), [&]( array::Array* a ) { size += a->footprint(); } );
    }
    return size;
}

// ----------------------------------------------------------------------------

void Table::dump( std::ostream& os ) const {
    values_.dump( os );
}

// ----------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
