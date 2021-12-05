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

#include "atlas/grid/detail/grid/CubedSphere.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/detail/grid/Unstructured.h"

namespace atlas {
namespace grid {

//----------------------------------------------------------------------------------------------------------------------

namespace {
void force_link() {
    static struct Link {
        Link() {
            GridFactoryBuilder<detail::grid::CubedSphere>();
            GridFactoryBuilder<detail::grid::Structured>();
            GridFactoryBuilder<detail::grid::Unstructured>();
        }
    } link;
}
}  // namespace

//----------------------------------------------------------------------------------------------------------------------

const GridImpl* GridFactory::build(const std::string& builder, const util::Config& config) {
    force_link();
    auto factory = get(builder);
    return factory->make(config);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
