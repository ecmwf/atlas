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

#include "eckit/memory/SharedPtr.h"

#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

namespace eckit {
class Configuration;
}

namespace atlas {
class FieldSet;
}

namespace atlas {

//------------------------------------------------------------------------------------------------------

class FunctionSpace {
public:
    using Implementation = functionspace::FunctionSpaceImpl;

private:
    eckit::SharedPtr<const Implementation> functionspace_;

public:
    FunctionSpace();
    FunctionSpace( const Implementation* );
    FunctionSpace( const FunctionSpace& );

    std::string type() const;
    operator bool() const;
    size_t footprint() const;
    std::string distribution() const;

    const Implementation* get() const { return functionspace_.get(); }

    atlas::Field createField( const eckit::Configuration& ) const;

    atlas::Field createField( const atlas::Field& ) const;
    atlas::Field createField( const atlas::Field&, const eckit::Configuration& ) const;

    template <typename DATATYPE>
    Field createField( const eckit::Configuration& ) const;

    template <typename DATATYPE>
    Field createField() const;

    void haloExchange( const FieldSet&, bool on_device = false ) const;
    void haloExchange( const Field&, bool on_device = false ) const;

    idx_t size() const { return functionspace_->size(); }
};


//------------------------------------------------------------------------------------------------------

}  // namespace atlas
