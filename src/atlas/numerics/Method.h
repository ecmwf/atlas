/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_numerics_Method_h
#define atlas_numerics_Method_h

#include <string>
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas {
namespace numerics {

/// @brief Method class
/// @note  Abstract base class
class Method : public eckit::Owned
{
public:
    Method() {}
    virtual ~Method() = 0;
    virtual std::string name() const = 0;
};

inline Method::~Method() {}

extern "C"
{
    void atlas__Method__delete (Method* This);
    const char* atlas__Method__name (Method* This);
}

} // namespace numerics
} // namespace atlas

#endif // atlas_numerics_Method_h
