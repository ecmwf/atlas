/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */


#pragma once

#include "atlas/interpolation/method/knn/GridBoxMethod.h"


namespace atlas {
namespace interpolation {
namespace method {


class GridBoxMaximum final : public GridBoxMethod {
public:
    using GridBoxMethod::GridBoxMethod;

private:
    void do_execute(const FieldSet& source, FieldSet& target, Metadata&) const override;
    void do_execute(const Field& source, Field& target, Metadata&) const override;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
