/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_FunctionSpace_h
#define atlas_FunctionSpace_h

#include <string>
#include "eckit/memory/Owned.h"

namespace atlas {
namespace functionspace {

/// @brief FunctionSpace class helps to interprete Fields.
/// @note  Abstract base class
class FunctionSpace : public eckit::Owned
{
public:
    FunctionSpace() {}
    virtual ~FunctionSpace() = 0;
    virtual std::string name() const = 0;
    eckit::SharedPtr<FunctionSpace const> shared_from_this() const;
    eckit::SharedPtr<FunctionSpace> shared_from_this();
    eckit::SharedPtr<FunctionSpace> ptr();
    eckit::SharedPtr<FunctionSpace const> ptr() const;
    eckit::SharedPtr<FunctionSpace const> cptr() const;
};

inline FunctionSpace::~FunctionSpace() {}

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C"
{
    void atlas__FunctionSpace__delete (FunctionSpace* This);
    const char* atlas__FunctionSpace__name (FunctionSpace* This);
}

//------------------------------------------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas

#endif // atlas_FunctionSpace_h
