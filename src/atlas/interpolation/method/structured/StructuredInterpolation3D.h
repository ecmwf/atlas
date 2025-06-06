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

#include "atlas/interpolation/method/Method.h"

#include <memory>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/grid/Vertical.h"

namespace atlas {
namespace interpolation {
namespace method {

/**
 * @class StructuredInterpolation3D
 *
 * Three-dimensional interpolation making use of Structure of grid.
 */

template <typename Kernel>
class StructuredInterpolation3D : public Method {
public:
    StructuredInterpolation3D(const Config& config);

    virtual ~StructuredInterpolation3D() override {}

    virtual void print(std::ostream&) const override;

protected:
    void setup(const FunctionSpace& source);

    virtual const FunctionSpace& source() const override { return source_; }

    virtual const FunctionSpace& target() const override { return target_; }

private:
    virtual void do_setup(const Grid& source, const Grid& target, const Cache&) override;

    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) override;

    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;

    virtual void do_setup(const FunctionSpace& source, const Field& target) override;

    virtual void do_setup(const FunctionSpace& source, const FieldSet& target) override;

    virtual void do_execute(const Field& src, Field& tgt, Metadata&) const override;

    virtual void do_execute(const FieldSet& src, FieldSet& tgt, Metadata&) const override;

    template <typename Value, int Rank>
    void execute_impl(const Kernel& kernel, const FieldSet& src, FieldSet& tgt) const;

    static double convert_units_multiplier(const Field& field);

protected:
    Field target_ghost_;
    Field target_lonlat_;
    Field target_vertical_;

    Field target_3d_;

    FieldSet target_xyz_;

    FunctionSpace source_;
    FunctionSpace target_;

    bool matrix_free_;
    bool limiter_;

    std::unique_ptr<Kernel> kernel_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas

#include "StructuredInterpolation3D.tcc"
