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

#include <iosfwd>
#include <string>

#include "atlas/util/ObjectHandle.h"

namespace eckit {
class Hash;
class PathName;
}  // namespace eckit
namespace atlas {
namespace util {
class Config;
}
}  // namespace atlas

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------

/// @class Spec
/// @brief Grid definition as configuration
struct Spec : DOXYGEN_HIDE( public util::ObjectHandle<util::Config> ) {
    using Handle::Handle;
    Spec() = default;
    Spec( const std::string& id );

    bool operator==( const Spec& other ) const { return uid() == other.uid(); }
    bool operator!=( const Spec& other ) const { return uid() != other.uid(); }

    std::string name() const;
    std::string type() const;
    std::string uid() const;

    /// Adds to the hash the information that makes this Spec unique
    void hash( eckit::Hash& h ) const;
};

//---------------------------------------------------------------------------------------------------------------------

struct SpecFactory {
    explicit SpecFactory( const std::string& id, const Spec::Implementation& );
    SpecFactory( const std::string& id, const eckit::PathName& );
    virtual ~SpecFactory() = default;

    static Spec::Implementation* create( const std::string& id );
    static void list( std::ostream& );

private:
    SpecFactory( const SpecFactory& ) = delete;
    SpecFactory& operator=( const SpecFactory& ) = delete;
    std::string name_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
