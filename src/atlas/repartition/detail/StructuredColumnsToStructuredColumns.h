#pragma once

#include "atlas/repartition/detail/RepartitionImpl.h"

namespace atlas {

  class Field;
  class FieldSet;
  class FunctionSpace;

  namespace functionspace {
    class StructuredColumns;
  }
}

namespace atlas {
  namespace repartition {

    using functionspace::StructuredColumns;

    // Concrete implementation of StructuredColumns to StructuredColumns.
    class StructuredColumnsToStructuredColumns : public RepartitionImpl {

    public:

      StructuredColumnsToStructuredColumns(
        const FunctionSpace& sourceFunctionSpace,
        const FunctionSpace& targetFunctionSpace);

      void execute(const Field &source, Field &target) const override;
      void execute(const FieldSet &source, FieldSet &target) const override;

    private:

      // FunctionSpaces recast to StructuredColumns.
      StructuredColumns& sourceStructuredColumns_;
      StructuredColumns& targetStructuredColumns_;

    };
  }
}
