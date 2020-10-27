/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/Interpolation.h"
namespace atlas {
namespace interpolation {

class MatrixCacheEntryOwned : public MatrixCacheEntry {
public:
    MatrixCacheEntryOwned( Matrix&& matrix ) : MatrixCacheEntry( &matrix_ ) {
        const_cast<Matrix&>( matrix_ ).swap( reinterpret_cast<Matrix&>( matrix ) );
    }

private:
    const Matrix matrix_;
};

class MatrixCacheEntryEmpty : public MatrixCacheEntryOwned {
public:
    MatrixCacheEntryEmpty() : MatrixCacheEntryOwned( Matrix{} ) {}
};

class MatrixCacheEntryShared : public MatrixCacheEntry {
public:
    MatrixCacheEntryShared( std::shared_ptr<const Matrix> matrix ) :
        MatrixCacheEntry( matrix.get() ), matrix_{matrix} {}

private:
    const std::shared_ptr<const Matrix> matrix_;
};


MatrixCache::MatrixCache( const Cache& c ) :
    Cache( c ), matrix_{dynamic_cast<const MatrixCacheEntry*>( c.get( MatrixCacheEntry::static_type() ) )} {}

MatrixCache::MatrixCache( Matrix&& m ) : MatrixCache( std::make_shared<MatrixCacheEntryOwned>( std::move( m ) ) ) {}

MatrixCache::MatrixCache( std::shared_ptr<const Matrix> m ) :
    MatrixCache( std::make_shared<MatrixCacheEntryShared>( m ) ) {}

MatrixCache::MatrixCache( const Matrix* m ) : MatrixCache( std::make_shared<MatrixCacheEntry>( m ) ) {}

MatrixCache::MatrixCache( const Interpolation& interpolation ) : MatrixCache( interpolation.get()->matrix_shared_ ) {}

MatrixCache::operator bool() const {
    return matrix_ && !matrix().empty();
}

const MatrixCache::Matrix& MatrixCache::matrix() const {
    ATLAS_ASSERT( matrix_ );
    return matrix_->matrix();
}

MatrixCache::MatrixCache( std::shared_ptr<InterpolationCacheEntry> entry ) :
    Cache( entry ), matrix_{dynamic_cast<const MatrixCacheEntry*>( entry.get() )} {}


}  // namespace interpolation
}  // namespace atlas
