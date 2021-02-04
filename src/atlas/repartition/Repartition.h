#pragma once

#include "atlas/repartition/detail/RepartitionImpl.h"
#include "atlas/util/ObjectHandle.h"

namespace atlas {
  class Field;
  class FieldSet;
  class FunctionSpace;
}

namespace atlas {

  // Object handle class for reparitioner
  class Repartition :
    public util::ObjectHandle<repartition::RepartitionImpl> {

  public:

    using Handle::Handle;

    Repartition () = default;

    // Set up repartitioner form source to target function space.
    Repartition (
      const FunctionSpace& sourceFunctionSpace,
      const FunctionSpace& targetFunctionSpace);

    void execute (const Field& sourceField, Field& targetField) const;
    void execute (
      const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const;

    FunctionSpace& getSourceFunctionSpace ();
    const FunctionSpace& getSourceFunctionSpace () const;
    FunctionSpace& getTargetFunctionSpace ();
    const FunctionSpace& getTargetFunctionSpace () const;

  };

}
