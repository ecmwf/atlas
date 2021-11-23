/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/redistribution/detail/RedistributionImpl.h"

namespace atlas {
namespace redistribution {
namespace detail {

void RedistributionImpl::setup( const FunctionSpace& source, const FunctionSpace& target ) {
    source_ = source;
    target_ = target;
    do_setup();
}

const FunctionSpace& RedistributionImpl::source() const {
    return source_;
}

const FunctionSpace& RedistributionImpl::target() const {
    return target_;
}

}  // namespace detail
}  // namespace redistribution
}  // namespace atlas
