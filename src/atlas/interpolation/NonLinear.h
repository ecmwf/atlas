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

#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {
class NonLinear;
}
}  // namespace interpolation
namespace util {
class Config;
}
}  // namespace atlas


namespace atlas {
namespace interpolation {


/**
 * @brief NonLinear class applies non-linear corrections to an interpolation matrix, given a field with missing values.
 * The interpolatation are re-weighted to factor those values out of the resulting field.
 */
struct NonLinear : DOXYGEN_HIDE( public util::ObjectHandle<nonlinear::NonLinear> ) {
    using Spec = util::Config;

    /// ctor
    using Handle::Handle;
    NonLinear();
    NonLinear( const util::Config& );
    NonLinear( const std::string& type, const util::Config& );

    /// bool operator
    using Handle::operator bool;  // (ensure this exists)
};


}  // namespace interpolation
}  // namespace atlas
