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

InterpolationCacheEntry::~InterpolationCacheEntry() = default;

Cache::Cache(std::shared_ptr<InterpolationCacheEntry> cache) {
    auto& new_cache = cache_[cache->type()];
    new_cache       = cache;
}

Cache::Cache(const Cache& other) {
    add(other);
}

Cache::Cache(const Interpolation& interpolation): Cache(interpolation.createCache()) {}

Cache::~Cache() = default;

void Cache::add(const Cache& other) {
    for (auto& entry : other.cache_) {
        auto& new_cache = cache_[entry.first];
        new_cache       = entry.second;
    }
}

MatrixCacheEntry::~MatrixCacheEntry() = default;

class MatrixCacheEntryOwned : public MatrixCacheEntry {
public:
    MatrixCacheEntryOwned(Matrix&& matrix): MatrixCacheEntry(&matrix_) {
        const_cast<Matrix&>(matrix_).swap(reinterpret_cast<Matrix&>(matrix));
    }

private:
    const Matrix matrix_;
};

class MatrixCacheEntryShared : public MatrixCacheEntry {
public:
    MatrixCacheEntryShared(std::shared_ptr<const Matrix> matrix): MatrixCacheEntry(matrix.get()), matrix_{matrix} {}

private:
    const std::shared_ptr<const Matrix> matrix_;
};


MatrixCache::MatrixCache(const Cache& c):
    Cache(c), matrix_{dynamic_cast<const MatrixCacheEntry*>(c.get(MatrixCacheEntry::static_type()))} {}

MatrixCache::MatrixCache(Matrix&& m): MatrixCache(std::make_shared<MatrixCacheEntryOwned>(std::move(m))) {}

MatrixCache::MatrixCache(std::shared_ptr<const Matrix> m): MatrixCache(std::make_shared<MatrixCacheEntryShared>(m)) {}

MatrixCache::MatrixCache(const Matrix* m): MatrixCache(std::make_shared<MatrixCacheEntry>(m)) {}

MatrixCache::MatrixCache(const Interpolation& interpolation): MatrixCache(Cache(interpolation)) {}

MatrixCache::operator bool() const {
    return matrix_ && !matrix().empty();
}

const MatrixCache::Matrix& MatrixCache::matrix() const {
    ATLAS_ASSERT(matrix_);
    return matrix_->matrix();
}

size_t MatrixCache::footprint() const {
    if (matrix_) {
        return matrix_->footprint();
    }
    return 0;
}

MatrixCache::MatrixCache(std::shared_ptr<InterpolationCacheEntry> entry):
    Cache(entry), matrix_{dynamic_cast<const MatrixCacheEntry*>(entry.get())} {}


IndexKDTreeCacheEntry::~IndexKDTreeCacheEntry() = default;

const IndexKDTreeCacheEntry::IndexKDTree& IndexKDTreeCacheEntry::tree() const {
    ATLAS_ASSERT(tree_);
    return tree_;
}

IndexKDTreeCache::IndexKDTreeCache(const Cache& c):
    Cache(c), tree_{dynamic_cast<const IndexKDTreeCacheEntry*>(c.get(IndexKDTreeCacheEntry::static_type()))} {}

IndexKDTreeCache::IndexKDTreeCache(const IndexKDTree& tree): Cache(std::make_shared<IndexKDTreeCacheEntry>(tree)) {
    tree_ = dynamic_cast<const IndexKDTreeCacheEntry*>(get(IndexKDTreeCacheEntry::static_type()));
}

IndexKDTreeCache::IndexKDTreeCache(const Interpolation& interpolation): IndexKDTreeCache(Cache(interpolation)) {}

IndexKDTreeCache::operator bool() const {
    return tree_;  //&& !tree().empty();
}

const IndexKDTreeCache::IndexKDTree& IndexKDTreeCache::tree() const {
    ATLAS_ASSERT(tree_);
    return tree_->tree();
}

size_t IndexKDTreeCache::footprint() const {
    if (tree_) {
        return tree_->footprint();
    }
    return 0;
}


}  // namespace interpolation
}  // namespace atlas
