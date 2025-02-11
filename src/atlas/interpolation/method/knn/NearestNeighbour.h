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

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/method/knn/KNearestNeighboursBase.h"

namespace atlas {
namespace interpolation {
namespace method {

class NearestNeighbour : public KNearestNeighboursBase {
public:
    NearestNeighbour(const Config& config): KNearestNeighboursBase(config) {}
    virtual ~NearestNeighbour() override {}

    virtual void print(std::ostream&) const override {}

protected:
    virtual const FunctionSpace& source() const override { return source_; }
    virtual const FunctionSpace& target() const override { return target_; }

private:
    using KNearestNeighboursBase::do_setup;

    /**
   * @brief Create an interpolant sparse matrix relating two (pre-partitioned)
   * meshes,
   * using nearest neighbour method
   * @param source functionspace containing source elements
   * @param target functionspace containing target points
   */
    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
    virtual void do_setup(const Grid& source, const Grid& target, const Cache&) override;
    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) override;

    FunctionSpace source_;
    FunctionSpace target_;
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
