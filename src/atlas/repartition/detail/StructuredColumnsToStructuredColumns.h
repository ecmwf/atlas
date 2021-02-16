#pragma once

#include <vector>

#include "atlas/functionspace/StructuredColumns.h"
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

    using idxPair = std::pair<idx_t, idx_t>;
    using idxPairVector = std::vector<idxPair>;

    // Helper struct for function space intersections.
    struct indexRange {
      idxPair jBeginEnd{};
      idxPairVector iBeginEnd{};
    };

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
      StructuredColumns sourceStructuredColumns_{};
      StructuredColumns targetStructuredColumns_{};

      // Function space intersection ranges.
      std::vector<indexRange> sendIntersections_{};
      std::vector<indexRange> receiveIntersections_{};

      // MPI count and displacement arrays for send and receive buffers.
      std::vector<int> sendCounts_{};
      std::vector<int> sendDisplacements_{};
      std::vector<int> receiveCounts_{};
      std::vector<int> receiveDisplacements_{};


    };
  }
}
