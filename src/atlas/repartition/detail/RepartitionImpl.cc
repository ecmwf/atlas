#include "atlas/parallel/mpi/mpi.h"
#include "atlas/repartition/detail/RepartitionImpl.h"

namespace atlas {
  namespace repartition {

    // Constructors/destructors.

    RepartitionImpl::RepartitionImpl (
      const FunctionSpace& source, const FunctionSpace& target) :
      sourceFunctionSpace_(source), targetFunctionSpace_(target) {}

    RepartitionImpl::~RepartitionImpl () {}

    // Getters.

    FunctionSpace& RepartitionImpl::getSourceFunctionSpace () {
      return sourceFunctionSpace_;
    }

    const FunctionSpace& RepartitionImpl::getSourceFunctionSpace () const {
      return sourceFunctionSpace_;
    }

    FunctionSpace& RepartitionImpl::getTargetFunctionSpace () {
      return targetFunctionSpace_;
    }

    const FunctionSpace& RepartitionImpl::getTargetFunctionSpace () const {
      return targetFunctionSpace_;
    }

  }
}
