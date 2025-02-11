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

#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/method/PointIndex3.h"
#include "atlas/mesh/Elements.h"

namespace atlas {
namespace interpolation {
namespace method {

class FiniteElement : public Method {
public:
    FiniteElement(const Config& config): Method(config) {
        config.get("max_fraction_elems_to_try", max_fraction_elems_to_try_);
        config.get("treat_failure_as_missing_value", treat_failure_as_missing_value_);
    }

    virtual ~FiniteElement() override {}

    virtual void print(std::ostream&) const override;

protected:
    /**
   * @brief Create an interpolant sparse matrix relating two (pre-partitioned)
   * meshes,
   * using elements as per the Finite Element Method and ray-tracing to
   * calculate
   * source mesh elements intersections (and interpolation weights) with target
   * grid
   * node-containing rays
   * @param meshSource mesh containing source elements
   * @param meshTarget mesh containing target points
   */
    void setup(const FunctionSpace& source);

    /**
   * Find in which element the point is contained by projecting (ray-tracing)
   * the
   * point to the nearest element(s), returning the (normalized) interpolation
   * weights
   */
    Triplets projectPointToElements(size_t ip, const ElemIndex3::NodeList& elems, std::ostream& failures_log) const;

    virtual const FunctionSpace& source() const override { return source_; }
    virtual const FunctionSpace& target() const override { return target_; }

private:
    using Method::do_setup;
    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
    virtual void do_setup(const Grid& source, const Grid& target, const Cache&) override;
    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) override;

protected:
    mesh::MultiBlockConnectivity* connectivity_;
    std::unique_ptr<array::ArrayView<double, 2>> icoords_;
    std::unique_ptr<array::ArrayView<double, 2>> ocoords_;
    std::unique_ptr<array::ArrayView<gidx_t, 1>> igidx_;

    Field target_lonlat_;
    Field target_xyz_;
    Field target_ghost_;

    FunctionSpace source_;
    FunctionSpace target_;

    bool treat_failure_as_missing_value_{true};
    double max_fraction_elems_to_try_{0.2};
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
