#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/repartition/detail/RepartitionImplFactory.h"
#include "atlas/repartition/detail/StructuredColumnsToStructuredColumns.h"

namespace atlas {
  namespace repartition {

    RepartitionImpl* RepartitionImplFactory::build (
      const FunctionSpace& sourceFunctionSpace,
      const FunctionSpace& targetFunctionSpace) {

      if (sourceFunctionSpace->type() == StructuredColumns::type() &&
          targetFunctionSpace->type() == StructuredColumns::type()) {

        std::cout << "Making StructuredColumnsToStructuredColumns" << std::endl;
        return new StructuredColumnsToStructuredColumns (
          sourceFunctionSpace, targetFunctionSpace);
      }
      else {

        std::cout << "Unknown source " + sourceFunctionSpace->type() << std::endl;
        std::cout << "Unknown target " + targetFunctionSpace->type() << std::endl;
        return nullptr;

      }
    }

  }
}

