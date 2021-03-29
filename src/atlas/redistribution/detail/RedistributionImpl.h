/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/util/Object.h"
#include "atlas/functionspace.h"

namespace atlas {
  class Field;
  class FieldSet;
  class FunctionSpace;
}

namespace atlas {
  namespace redistribution {
    namespace detail {

      /// \brief  Abstract base class for redistributor implementation.
      class RedistributionImpl : public util::Object {

      public:

        /// \brief  Virtual destructor.
        virtual ~RedistributionImpl() = 0;

        /// \brief  Maps source field to target field.
        virtual void execute(
          const Field& sourceField, Field& targetField) const = 0;

        /// \brief  Maps source field set to target field set.
        virtual void execute(
          const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const = 0;

        /// \brief  Get reference to source function space.
        FunctionSpace& source();

        /// \brief  Get const reference to source function space.
        const FunctionSpace& source() const;

        /// \brief  Get reference to target function space.
        FunctionSpace& target();

        /// \brief  Get const reference to target function space.
        const FunctionSpace& target() const;

      protected:

        RedistributionImpl (
          const FunctionSpace& source, const FunctionSpace& target);

      private:

        FunctionSpace sourceFunctionSpace_;
        FunctionSpace targetFunctionSpace_;

      };

    }
  }
}
