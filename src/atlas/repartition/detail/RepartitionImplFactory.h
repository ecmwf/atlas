/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace atlas {

  class FunctionSpace;

  namespace repartition {
    class RepartitionImpl;
  }
}

namespace atlas {
  namespace repartition{

    /// \brief  Factory class to select correct concrete repartitioner.
    class RepartitionImplFactory {

    public:

      /// \brief  Selection based on source and target function spaces.
      static RepartitionImpl* build(
        const FunctionSpace& sourceFunctionSpace,
        const FunctionSpace& targetFunctionSpace);
    };

  }
}
