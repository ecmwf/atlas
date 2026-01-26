/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <map>
#include <memory>
#include <string>

#include "eckit/filesystem/PathName.h"
#include "eckit/io/Buffer.h"


#include "atlas/linalg/sparse.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/KDTree.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Interpolation;
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace interpolation {

//-----------------------------------------------------------------------------

class InterpolationCacheEntry {
public:
    virtual ~InterpolationCacheEntry();
    virtual size_t footprint() const = 0;
    virtual std::string type() const = 0;
};

//-----------------------------------------------------------------------------

class Cache {
public:
    Cache() = default;
    Cache(const Cache& other);
    Cache(const Cache& other, const std::string& filter);
    Cache(const Interpolation&);
    Cache& operator=(const Cache& other);
    Cache& operator=(Cache&& other);

    operator bool() const { return not cache_.empty(); }
    virtual ~Cache();
    size_t footprint() const {
        size_t footprint{0};
        for (auto& entry : cache_) {
            footprint += entry.second->footprint();
        }
        return footprint;
    }
    void add(const Cache&);

protected:
    Cache(std::shared_ptr<InterpolationCacheEntry> cache);

public:
    const InterpolationCacheEntry* get(const std::string& type) const {
        auto it = cache_.find(type);
        if (it != cache_.end()) {
            return it->second.get();
        }
        return nullptr;
    }

private:
    std::map<std::string, std::shared_ptr<InterpolationCacheEntry>> cache_;
};

//-----------------------------------------------------------------------------

class MatrixCacheEntry : public InterpolationCacheEntry {
public:
    using Matrix = atlas::linalg::SparseMatrixStorage;
    ~MatrixCacheEntry() override;
    MatrixCacheEntry(const Matrix* matrix, const std::string& uid = ""): matrix_{matrix}, uid_(uid) {
        ATLAS_ASSERT(matrix_ != nullptr);
    }
    const Matrix& matrix() const { return *matrix_; }
    const std::string& uid() const { return uid_; }
    size_t footprint() const override { return matrix_->footprint(); }
    operator bool() const { return not matrix_->empty(); }
    static std::string static_type() { return "Matrix"; }
    std::string type() const override { return static_type(); }

private:
    const Matrix* matrix_;
    const std::string uid_;
};

//-----------------------------------------------------------------------------

class MatrixCache final : public Cache {
public:
    using Matrix = MatrixCacheEntry::Matrix;

public:
    MatrixCache() = default;
    MatrixCache(const Cache& c);
    MatrixCache(Matrix&& m);
    MatrixCache(std::shared_ptr<const Matrix> m, const std::string& uid = "");
    MatrixCache(const Matrix* m);
    MatrixCache(const Interpolation&);

    operator bool() const;
    const Matrix& matrix() const;
    const std::string& uid() const;
    size_t footprint() const;

private:
    MatrixCache(std::shared_ptr<InterpolationCacheEntry> entry);
    const MatrixCacheEntry* matrix_{nullptr};
};

//----------------------------------------------------------------------------------------------------------------------

class IndexKDTreeCacheEntry : public InterpolationCacheEntry {
public:
    using IndexKDTree = util::IndexKDTree;
    IndexKDTreeCacheEntry(const IndexKDTree& tree): tree_{tree} { ATLAS_ASSERT(tree_); }
    virtual ~IndexKDTreeCacheEntry() override;
    const IndexKDTree& tree() const;
    size_t footprint() const override { return tree().footprint(); }
    operator bool() const { return bool(tree_); }
    static std::string static_type() { return "IndexKDTree"; }
    std::string type() const override { return static_type(); }

private:
    const IndexKDTree tree_;
};

class IndexKDTreeCache final : public Cache {
public:
    using IndexKDTree = IndexKDTreeCacheEntry::IndexKDTree;

public:
    IndexKDTreeCache() = default;
    IndexKDTreeCache(const Cache& c);
    IndexKDTreeCache(const IndexKDTreeCache& c) : IndexKDTreeCache(Cache(c)) {}
    IndexKDTreeCache(const IndexKDTree&);
    IndexKDTreeCache(const Interpolation&);
    operator bool() const;
    const IndexKDTree& tree() const;
    size_t footprint() const;

private:
    const IndexKDTreeCacheEntry* tree_{nullptr};
};

}  // namespace interpolation
}  // namespace atlas
