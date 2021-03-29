/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/redistribution/detail/RedistributionImpl.h"
#include "atlas/util/ObjectHandle.h"

namespace atlas {
  class Field;
  class FieldSet;
  class FunctionSpace;
}

namespace atlas {

  /// \brief    Base redistributer class.
  ///
  /// \details  Class to map two function spaces with the same grid but
  ///           different partitioners.
  class Redistribution :
    public util::ObjectHandle<redistribution::detail::RedistributionImpl> {

  public:

    using Handle::Handle;

    /// \brief    Constructs and initialises the redistributor.
    ///
    /// \details  Initialises class to copy fields from a source function space
    ///           to fields from a target functionspace. The grids of source and
    ///           target function space must match.
    ///
    /// \param[in]  sourceFunctionSpace  Function space of source fields.
    /// \param[in]  targetFunctionSpace  Function space of target fields.
    Redistribution(
      const FunctionSpace& sourceFunctionSpace,
      const FunctionSpace& targetFunctionSpace);

    /// \brief    Maps source field to target field.
    ///
    /// \details  Transfers data from source field to target field. Function
    ///           space of source field must match sourceFunctionSpace. Same
    ///           applies to target field.
    ///
    /// \param[in]  sourceField  input field matching sourceFunctionSpace.
    /// \param[out] targetField  output field matching targetFunctionSpace.
    void execute(const Field& sourceField, Field& targetField) const;

    /// \brief    Redistributes source field set to target fields set.
    ///
    /// \details  Transfers source field set to target field set via multiple
    ///           invocations of execute(sourceField, targetField).
    ///
    /// \param[in]  sourceFieldSet  input field set.
    /// \param[out] targetFieldSet  output field set.
    void execute(
      const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const;

    /// \brief  Get reference to source function space.
    FunctionSpace& source();

    /// \brief  Get const reference to source function space.
    const FunctionSpace& source() const;

    /// \brief  Get reference to taget function space.
    FunctionSpace& target();

    /// \brief  Get const reference to target function space.
    const FunctionSpace& target() const;

  };

}
