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

#include <string>

#include "atlas/grid/detail/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

class FocusSpacing : public Spacing {
public:
    // constructor
    FocusSpacing(const eckit::Parametrisation& p);

    // class name
    static std::string static_type() { return "focus"; }
    virtual std::string type() const { return static_type(); }

    virtual Spec spec() const;

private:
    double focus_factor_;
    double start_;
    double end_;
};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
