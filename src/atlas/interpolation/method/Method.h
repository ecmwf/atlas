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

#include <iosfwd>
#include <string>
#include <vector>

#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/NonLinear.h"
#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Object.h"
#include "eckit/config/Configuration.h"
#include "eckit/linalg/SparseMatrix.h"
#include "atlas/linalg/sparse/MakeSparseMatrixStorageEckit.h"

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
class Grid;
}  // namespace atlas

namespace atlas {
namespace test {
class Access;
}
}  // namespace atlas

namespace atlas {
namespace interpolation {

class Method : public util::Object {
public:
    using Config   = eckit::Parametrisation;
    using Metadata = util::Metadata;

    Method(const Config&);
    virtual ~Method() {}

    /**
     * @brief Setup the interpolator relating two functionspaces
     * @param source functionspace containing source elements
     * @param target functionspace containing target points
     */
    void setup(const FunctionSpace& source, const FunctionSpace& target);
    void setup(const Grid& source, const Grid& target);
    void setup(const FunctionSpace& source, const Field& target);
    void setup(const FunctionSpace& source, const FieldSet& target);
    void setup(const Grid& source, const Grid& target, const Cache&);
    void setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&);

    Metadata execute(const FieldSet& source, FieldSet& target) const;
    Metadata execute(const Field& source, Field& target) const;

    /**
     * @brief execute_adjoint
     * @param source - it is either a FieldSet or a Field
     * @param target - it is either a FieldSet or a Field
     *                 Note that formally in an adjoint operation of this
     *                 type we should be setting the values in the target
     *                 to zero. This is not done for efficiency reasons and
     *                 because in most cases it is not necessary.
     */
    Metadata execute_adjoint(FieldSet& source, const FieldSet& target) const;
    Metadata execute_adjoint(Field& source, const Field& target) const;

    virtual void print(std::ostream&) const = 0;

    virtual const FunctionSpace& source() const = 0;
    virtual const FunctionSpace& target() const = 0;

    virtual interpolation::Cache createCache() const;

protected:
    virtual void do_execute(const FieldSet& source, FieldSet& target, Metadata&) const;
    virtual void do_execute(const Field& source, Field& target, Metadata&) const;

    virtual void do_execute_adjoint(FieldSet& source, const FieldSet& target, Metadata&) const;
    virtual void do_execute_adjoint(Field& source, const Field& target, Metadata&) const;

    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;
    using Matrix   = atlas::linalg::SparseMatrixStorage;

    static void normalise(Triplets& triplets);

    void haloExchange(const FieldSet&) const;
    void haloExchange(const Field&) const;

    void adjointHaloExchange(const FieldSet&) const;
    void adjointHaloExchange(const Field&) const;

    // NOTE : Matrix-free or non-linear interpolation operators do not have matrices, so do not expose here
    friend class atlas::test::Access;
    friend class interpolation::MatrixCache;

protected:

    void setMatrix(Matrix&& m, const std::string& uid = "") {
        if (not matrix_shared_) {
            matrix_shared_ = std::make_shared<Matrix>();
        }
        *matrix_shared_ = std::move(m);
        matrix_cache_ = interpolation::MatrixCache(matrix_shared_, uid);
        matrix_       = &matrix_cache_.matrix();

    }

    void setMatrix(interpolation::MatrixCache matrix_cache) {
        ATLAS_ASSERT(matrix_cache);
        matrix_cache_ = matrix_cache;
        matrix_       = &matrix_cache_.matrix();
        matrix_shared_.reset();
    }

    void setMatrix(eckit::linalg::SparseMatrix&& m, const std::string& uid = "") {
        setMatrix( linalg::make_sparse_matrix_storage(std::move(m)), uid );
    }

    void setMatrix(std::size_t rows, std::size_t cols, const Triplets& triplets, const std::string& uid = "") {
        setMatrix( eckit::linalg::SparseMatrix{rows, cols, triplets}, uid);
    }

    bool matrixAllocated() const { return matrix_shared_.use_count(); }

    const Matrix& matrix() const { return *matrix_; }

    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target) = 0;
    virtual void do_setup(const Grid& source, const Grid& target, const Cache&)     = 0;
    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) = 0;
    virtual void do_setup(const FunctionSpace& source, const Field& target);
    virtual void do_setup(const FunctionSpace& source, const FieldSet& target);

    void check_compatibility(const Field& src, const Field& tgt, const Matrix& W) const;

private:
    template <typename Value>
    void interpolate_field(const Field& src, Field& tgt, const Matrix&) const;

    template <typename Value>
    void interpolate_field_rank1(const Field& src, Field& tgt, const Matrix&) const;

    template <typename Value>
    void interpolate_field_rank2(const Field& src, Field& tgt, const Matrix&) const;

    template <typename Value>
    void interpolate_field_rank3(const Field& src, Field& tgt, const Matrix&) const;

    template <typename Value>
    void adjoint_interpolate_field(Field& src, const Field& tgt, const Matrix&) const;

    template <typename Value>
    void adjoint_interpolate_field_rank1(Field& src, const Field& tgt, const Matrix&) const;

    template <typename Value>
    void adjoint_interpolate_field_rank2(Field& src, const Field& tgt, const Matrix&) const;

    template <typename Value>
    void adjoint_interpolate_field_rank3(Field& src, const Field& tgt, const Matrix&) const;

private:
    const Matrix* matrix_ = nullptr;
    std::shared_ptr<Matrix> matrix_shared_;
    interpolation::MatrixCache matrix_cache_;
    NonLinear nonLinear_;
    std::string linalg_backend_;
    Matrix matrix_transpose_;

protected:
    bool adjoint_{false};
    bool allow_halo_exchange_{true};
    std::vector<idx_t> missing_;
};

}  // namespace interpolation
}  // namespace atlas
