#include "atlas/field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/repartition/detail/StructuredColumnsToStructuredColumns.h"

namespace atlas {
  namespace repartition {

    StructuredColumnsToStructuredColumns::StructuredColumnsToStructuredColumns (
      const FunctionSpace& sourceFunctionSpace,
      const FunctionSpace& targetFunctionSpace) :
      RepartitionImpl (sourceFunctionSpace, targetFunctionSpace),
      sourceStructuredColumns_ (
        dynamic_cast<StructuredColumns&> (getSourceFunctionSpace ())),
      targetStructuredColumns_ (
        dynamic_cast<StructuredColumns&> (getTargetFunctionSpace ())) {

      std::cout << "Constructor on MPI rank " << mpi::rank () << std::endl;

      return;
    }

    void StructuredColumnsToStructuredColumns::execute (
      const Field& sourceField, Field& targetField) const {

      std::cout << "Field execute on MPI rank " << mpi::rank () << std::endl;

      return;
    }

    void StructuredColumnsToStructuredColumns::execute (
      const FieldSet& sourceField, FieldSet& targetField) const {

      std::cout << "FieldSet execute on MPI rank " << mpi::rank () << std::endl;

      return;
    }

  }
}
