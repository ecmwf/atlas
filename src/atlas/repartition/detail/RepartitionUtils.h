/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <typeinfo>

#include "atlas/field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

#include "eckit/exception/Exceptions.h"

namespace atlas {
  namespace repartition {

    using functionspace::FunctionSpaceImpl;

    /// \brief  Check function space can be cast to FunctionSpaceType.
    template <typename FunctionSpaceType>
    void tryCast(const FunctionSpaceImpl* const functionSpacePtr,
      const std::string& varName) {

      // Check if cast failed.
      if (!(functionSpacePtr->cast<FunctionSpaceType>())) {
         throw eckit::BadCast("Cannot cast " + varName + " to " +
         typeid(FunctionSpaceType).name() , Here());
      }

      return;
    }

    /// \brief  Check grids associated with two function spaces match.
    template <typename FunctionSpaceType>
    void checkGrids(const FunctionSpaceImpl* const functionSpacePtrA,
      const FunctionSpaceImpl* const functionSpacePtrB,
      const std::string& varNameA, const std::string& varNameB) {

      // Cast function spaces.
      const auto* const castPtrA = functionSpacePtrA->cast<FunctionSpaceType>();
      const auto* const castPtrB = functionSpacePtrB->cast<FunctionSpaceType>();

      // Check casts.
      tryCast<FunctionSpaceType>(castPtrA, varNameA);
      tryCast<FunctionSpaceType>(castPtrB, varNameB);

      // Check grids match.
      const auto gridNameA = castPtrA->grid().name();
      const auto gridNameB = castPtrB->grid().name();

      if (gridNameA != gridNameB) {

        throw eckit::BadValue("Grids do not match.\n" +
          varNameA + "->grid() is type " + gridNameA + "\n" +
          varNameB + "->grid() is type " + gridNameB, Here());
      }

      return;
    }

    /// \brief  Check fields have the same data type.
    void checkFieldDataType(const Field& fieldA, const Field& fieldB,
      const std::string& varnameA, const std::string& varnameB) {

      // Check fields have the same data type.
      if (fieldA.datatype() != fieldB.datatype()) {
        throw eckit::BadValue("Fields have different data types.\n" +
          varnameA + " has data type " + fieldA.datatype().str() + "\n" +
          varnameB + " has data type " + fieldB.datatype().str(), Here());
      }

      return;
    }

    /// \brief  Check field sets the same size.
    void checkFieldSetSize(const FieldSet& fieldSetA, const FieldSet& fieldSetB,
      const std::string& varnameA, const std::string& varnameB) {

      // Check fields have the same data type.
      if (fieldSetA.size() != fieldSetB.size()) {
        throw eckit::BadValue("Field sets have different data types.\n" +
          varnameA + " has size " + std::to_string(fieldSetA.size()) + "\n" +
          varnameB + " has size " + std::to_string(fieldSetB.size()), Here());
      }

      return;
    }
  }
}
