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

      // Constructors/destructors.

      RedistributionImpl::RedistributionImpl (
        const FunctionSpace& source, const FunctionSpace& target) :
        sourceFunctionSpace_(source), targetFunctionSpace_(target) {}

      RedistributionImpl::~RedistributionImpl() {}

      // Getters.

      FunctionSpace& RedistributionImpl::source() {
        return sourceFunctionSpace_;
      }

      const FunctionSpace& RedistributionImpl::source() const {
        return sourceFunctionSpace_;
      }

      FunctionSpace& RedistributionImpl::target() {
        return targetFunctionSpace_;
      }

      const FunctionSpace& RedistributionImpl::target() const {
        return targetFunctionSpace_;
      }

    }
  }
}
