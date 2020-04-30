/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "atlas/grid/detail/spacing/SpacingFactory.h"

namespace atlas {
namespace grid {
namespace spacing {

//----------------------------------------------------------------------------------------------------------------------


void force_link() {
    static struct Link { Link() = default; } link;
    []( const Link& ) {}( link );  // disable unused warnings
}

//----------------------------------------------------------------------------------------------------------------------

const Spacing* SpacingFactory::build( const std::string& builder ) {
    return build( builder, util::NoConfig() );
}

const Spacing* SpacingFactory::build( const std::string& builder, const eckit::Parametrisation& param ) {
    force_link();
    auto factory = get( builder );
    return factory->make( param );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
