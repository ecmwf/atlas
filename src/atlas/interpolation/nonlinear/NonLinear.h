/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include <string>
#include <type_traits>

#include "eckit/config/Parametrisation.h"
#include "eckit/linalg/SparseMatrix.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Factory.h"
#include "atlas/util/ObjectHandle.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


/**
 * @brief NonLinear class applies non-linear corrections to an interpolation matrix, given a field with missing values.
 * The interpolatation are re-weighted to factor those values out of the resulting field.
 */
class NonLinear : public util::Object {
public:
    using Config = eckit::Parametrisation;
    using Matrix = eckit::linalg::SparseMatrix;
    using Scalar = eckit::linalg::Scalar;
    using Size   = eckit::linalg::Size;

    /**
     * @brief ctor
     */
    NonLinear() = default;

    /**
     * @brief dtor
     */
    virtual ~NonLinear() = default;

    /**
     * @bried if NonLinear applies to field
     * @param [in] f field
     * @return if NonLinear applies to field
     */
    virtual bool applicable(const Field& f) const = 0;

    /**
     * @brief Apply non-linear corrections to interpolation matrix
     * @param [inout] W interpolation matrix
     * @param [in] f field with missing values information
     * @return if W was modified
     */
    virtual bool execute(Matrix& W, const Field& f) const = 0;

protected:
    template <typename Value, int Rank>
    static array::ArrayView<typename std::add_const<Value>::type, Rank> make_view_field_values(const Field& field) {
        ATLAS_ASSERT(field);
        ATLAS_ASSERT_MSG(
            field.datatype().kind() == array::DataType::kind<Value>(),
            "Field(name:" + field.name() + ",DataType:" + field.datatype().str() + ") is not of required DataType");
        return array::make_view<typename std::add_const<Value>::type, Rank>(field);
    }
};


class NonLinearFactory : public util::Factory<NonLinearFactory> {
public:
    using Config = NonLinear::Config;

    static std::string className() { return "NonLinearFactory"; }
    static const NonLinear* build(const std::string&, const Config&);
    using Factory::Factory;

private:
    virtual const NonLinear* make(const Config&) = 0;
};


template <class T>
class NonLinearFactoryBuilder : public NonLinearFactory {
private:
    virtual const NonLinear* make(const Config& /*config*/) override { return new T(/*config*/); }

public:
    using NonLinearFactory::NonLinearFactory;
};


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
