/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/redistribution/detail/RedistributionImplFactory.h"
#include "atlas/redistribution/detail/StructuredColumnsToStructuredColumns.h"

#include "eckit/exception/Exceptions.h"

namespace atlas {
  namespace redistribution {
    namespace detail {

      RedistributionImpl* RedistributionImplFactory::build (
        const FunctionSpace& sourceFunctionSpace,
        const FunctionSpace& targetFunctionSpace) {


        if (sourceFunctionSpace->type() ==
          functionspace::StructuredColumns::type() &&
          targetFunctionSpace->type() ==
          functionspace::StructuredColumns::type()) {

          return new StructuredColumnsToStructuredColumns (
            sourceFunctionSpace, targetFunctionSpace);
        }
        else {

          throw eckit::NotImplemented(
            "No implementation for source function space " +
            sourceFunctionSpace->type() + " and target functionspace " +
            targetFunctionSpace->type(), Here());
          return nullptr;

        }
      }
    }
  }
}

