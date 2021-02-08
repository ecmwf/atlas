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

      // Helper struct for function space intersections.
      struct intersection {
        idx_t jBegin;
        idx_t jEnd;
        std::vector<idx_t> iBegin;
        std::vector<idx_t> iEnd;
      };

      // FunctionSpaces recast to StructuredColumns.
      StructuredColumns* sourceStructuredColumns_;
      StructuredColumns* targetStructuredColumns_;


    };
  }
}
