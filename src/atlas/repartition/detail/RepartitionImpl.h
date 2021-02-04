#pragma once

#include "atlas/util/Object.h"
#include "atlas/functionspace.h"

namespace atlas {
  class Field;
  class FieldSet;
  class FunctionSpace;
}

namespace atlas {
  namespace repartition {

    // Abstract class for repartitioner implementation.
    class RepartitionImpl : public util::Object {

    public:

      virtual ~RepartitionImpl () = 0;

      virtual void execute (
        const Field& sourceField, Field& targetField) const = 0;
      virtual void execute (
        const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const = 0;

      FunctionSpace& getSourceFunctionSpace ();
      const FunctionSpace& getSourceFunctionSpace () const;
      FunctionSpace& getTargetFunctionSpace ();
      const FunctionSpace& getTargetFunctionSpace () const;

    protected:

      RepartitionImpl (
        const FunctionSpace& source, const FunctionSpace& target);

    private:

      FunctionSpace sourceFunctionSpace_;
      FunctionSpace targetFunctionSpace_;

    };

  }
}
