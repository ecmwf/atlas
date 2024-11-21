#pragma once

#include <memory>
#include <vector>

#include "atlas/array.h"

#include "eckit/linalg/SparseMatrix.h"
#include "eckit/linalg/Triplet.h"
#include "eckit/linalg/types.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

/// Sparse matrix in CRS (compressed row storage) format
class SparseMatrix {
public:
    using Scalar = eckit::linalg::Scalar;
    using Index = eckit::linalg::Index;
    using Size = eckit::linalg::Size;
    using iterator = eckit::linalg::SparseMatrix::iterator;
    using const_iterator = eckit::linalg::SparseMatrix::const_iterator;

public:
    // -- Constructors

    /// Default constructor, empty matrix
    SparseMatrix();

    /// Constructor from triplets
    SparseMatrix(Size nrows, Size ncols, const std::vector<eckit::linalg::Triplet>&);

    /// Copy constructor
    SparseMatrix(const SparseMatrix&);

    /// Assignment operator (allocates and copies data)
    SparseMatrix& operator=(const SparseMatrix&);

public:
    void swap(SparseMatrix&);

    /// @returns number of rows
    Size rows() const { return host_matrix_.rows(); }

    /// @returns number of columns
    Size cols() const { return host_matrix_.cols(); }

    /// @returns number of non-zeros
    Size nonZeros() const { return host_matrix_.nonZeros(); }

    /// @returns true if this matrix does not contain non-zero entries
    bool empty() const { return host_matrix_.empty(); }

    /// @returns footprint of the matrix in memory
    size_t footprint() const;

    /// Prune entries, in-place, with exactly the given value
    SparseMatrix& prune(Scalar = 0);

    /// Transpose matrix in-place
    SparseMatrix& transpose();

    SparseMatrix& setIdentity(Size nrows, Size ncols);

public:
    void updateDevice() const { 
        outer_->updateDevice();
        inner_->updateDevice();
        data_->updateDevice();
    }

    void updateHost() const { 
        outer_->updateHost();
        inner_->updateHost();
        data_->updateHost();
    }

    bool hostNeedsUpdate() const { 
        return outer_->hostNeedsUpdate() ||
               inner_->hostNeedsUpdate() ||
               data_->hostNeedsUpdate();
    }

    bool deviceNeedsUpdate() const { 
        return outer_->deviceNeedsUpdate() ||
               inner_->deviceNeedsUpdate() ||
               data_->deviceNeedsUpdate();
    }

    void setHostNeedsUpdate(bool v) const { 
        outer_->setHostNeedsUpdate(v);
        inner_->setHostNeedsUpdate(v);
        data_->setHostNeedsUpdate(v);
    }

    void setDeviceNeedsUpdate(bool v) const {
        outer_->setDeviceNeedsUpdate(v);
        inner_->setDeviceNeedsUpdate(v);
        data_->setDeviceNeedsUpdate(v);
    }

    bool deviceAllocated() const {
        return outer_->deviceAllocated() &&
               inner_->deviceAllocated() &&
               data_->deviceAllocated();
    }

    void allocateDevice() {
        outer_->allocateDevice();
        inner_->allocateDevice();
        data_->allocateDevice();
    }

    void deallocateDevice() { 
        outer_->deallocateDevice();
        inner_->deallocateDevice();
        data_->deallocateDevice();
    }

    const eckit::linalg::SparseMatrix& host_matrix() const { return host_matrix_; }

    // eckit::linalg::SparseMatrix& host_matrix() { return host_matrix_; }

    const Scalar* data() const { return host_matrix_.data(); }
    
    const Index* outer() const { return host_matrix_.outer(); }

    const Index* inner() const { return host_matrix_.inner(); }

    const Scalar* host_data() const { return host_matrix_.data(); }
    
    const Index* host_outer() const { return host_matrix_.outer(); }

    const Index* host_inner() const { return host_matrix_.inner(); }

    const Scalar* device_data() const { return data_->device_data<Scalar>(); }
    
    const Index* device_outer() const { return outer_->device_data<Index>(); }

    const Index* device_inner() const { return inner_->device_data<Index>(); }

public:  // iterators
    
    /// const iterators to begin/end of row
    const_iterator begin(Size row) const { return host_matrix_.begin(row); }
    const_iterator end(Size row) const { return host_matrix_.end(row); }

    /// const iterators to begin/end of matrix
    const_iterator begin() const { return host_matrix_.begin(); }
    const_iterator end() const { return host_matrix_.end(); }

    /// iterators to begin/end of row
    iterator begin(Size row) { return host_matrix_.begin(row); }
    iterator end(Size row) { return host_matrix_.end(row); }

    /// const iterators to begin/end of matrix
    iterator begin() { return host_matrix_.begin(); }
    iterator end() { return host_matrix_.end(); }

private:
    void resetDeviceStorage();

private:
    eckit::linalg::SparseMatrix host_matrix_;
    std::unique_ptr<atlas::array::Array> outer_;
    std::unique_ptr<atlas::array::Array> inner_;
    std::unique_ptr<atlas::array::Array> data_;
};

//----------------------------------------------------------------------------------------------------------------------
}  // namespace linalg
}  // namespace atlas