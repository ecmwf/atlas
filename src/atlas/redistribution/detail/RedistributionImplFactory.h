/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace atlas {

  class FunctionSpace;

  namespace redistribution {
    namespace detail {
      class RedistributionImpl;
    }
  }
}

namespace atlas {
  namespace redistribution{
    namespace detail {

      /// \brief  Factory class to select correct concrete redistributor.
      class RedistributionImplFactory {

      public:

        /// \brief  Selection based on source and target function spaces.
        static RedistributionImpl* build(
          const FunctionSpace& sourceFunctionSpace,
          const FunctionSpace& targetFunctionSpace);
      };
    }
  }
}
