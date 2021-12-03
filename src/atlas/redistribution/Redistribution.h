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
}  // namespace atlas

namespace atlas {

/// \brief    Base redistributer class.
///
/// \details  Class to map two function spaces with the same grid but
///           different partitioners.
class Redistribution : public util::ObjectHandle<redistribution::detail::RedistributionImpl> {
public:
    using Handle::Handle;

    /// \brief    Empty default constructor.
    Redistribution();

    /// \brief    Constructs and initialises the redistributor.
    ///
    /// \details  Initialises class to copy fields from a source function space
    ///           to fields from a target functionspace. The grids of source and
    ///           target function space must match.
    ///
    /// \param[in]  source  Function space of source fields.
    /// \param[in]  target  Function space of target fields.
    /// \param[in]  config  Config to specify "type" of redistribution.
    Redistribution(const FunctionSpace& source, const FunctionSpace& target,
                   const util::Config& config = util::Config());

    /// \brief    Maps source field to target field.
    ///
    /// \details  Transfers data from source field to target field. Function
    ///           space of source field must match sourceFunctionSpace. Same
    ///           applies to target field.
    ///
    /// \param[in]  source  input field matching sourceFunctionSpace.
    /// \param[out] target  output field matching targetFunctionSpace.
    void execute(const Field& source, Field& target) const;

    /// \brief    Redistributes source field set to target fields set.
    ///
    /// \details  Transfers source field set to target field set via multiple
    ///           invocations of execute(sourceField, targetField).
    ///
    /// \param[in]  source  input field set.
    /// \param[out] target  output field set.
    void execute(const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const;

    /// \brief  Get const reference to source function space.
    const FunctionSpace& source() const;

    /// \brief  Get const reference to target function space.
    const FunctionSpace& target() const;
};

}  // namespace atlas
