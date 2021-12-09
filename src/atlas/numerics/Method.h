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

#include "atlas/util/Object.h"

namespace atlas {
namespace numerics {

/// @brief Method class
/// @note  Abstract base class
class Method : public util::Object {
public:
    Method() {}
    virtual ~Method()                       = 0;
    virtual const std::string& name() const = 0;
};

inline Method::~Method() {}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern "C" {
void atlas__Method__delete(Method* This);
const char* atlas__Method__name(Method* This);
}
#endif

}  // namespace numerics
}  // namespace atlas
