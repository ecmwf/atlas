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

#include "atlas/interpolation/nonlinear/MissingValue.h"
#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"


namespace atlas {
namespace util {
class Config;
}
}  // namespace atlas


namespace atlas {
namespace interpolation {


/// @brief Missing values indicator object
struct MissingValue : DOXYGEN_HIDE( public util::ObjectHandle<nonlinear::MissingValue> ) {
    using Spec   = util::Config;
    using Config = nonlinear::MissingValue::Config;

    using Handle::Handle;
    MissingValue();
    MissingValue( const Config& );
    MissingValue( const std::string& type, const Config& );

    using Handle::operator bool;  // (ensure this exists)

    bool operator()( const double& ) const;

    const nonlinear::MissingValue& ref() const;
};


}  // namespace interpolation
}  // namespace atlas
