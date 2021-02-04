#pragma once

namespace atlas {

  class FunctionSpace;

  namespace repartition {
    class RepartitionImpl;
  }
}

namespace atlas {
  namespace repartition{

    // Factory class to select correct concrete repartitioner.
    class RepartitionImplFactory {

    public:

      // Selection based on source and target FunctionSpaces.
      static RepartitionImpl* build (
        const FunctionSpace& sourceFunctionSpace,
        const FunctionSpace& targetFunctionSpace);
    };

  }
}
