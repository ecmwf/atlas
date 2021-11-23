/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "atlas/redistribution/detail/RedistributeGeneric.h"
#include "atlas/redistribution/detail/RedistributeStructuredColumns.h"
#include "atlas/redistribution/detail/RedistributionImpl.h"
#include "atlas/redistribution/detail/RedistributionImplFactory.h"

namespace atlas {
namespace redistribution {
namespace detail {

//----------------------------------------------------------------------------------------------------------------------

void force_link() {
    static struct Link {
        Link() {
            RedistributionImplBuilder<redistribution::detail::RedistributeGeneric>();
            RedistributionImplBuilder<redistribution::detail::RedistributeStructuredColumns>();
        }
    } link;
}

//----------------------------------------------------------------------------------------------------------------------

redistribution::detail::RedistributionImpl* RedistributionImplFactory::build( const std::string& builder ) {
    force_link();
    auto factory = get( builder );
    return factory->make();
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace redistribution
}  // namespace atlas
