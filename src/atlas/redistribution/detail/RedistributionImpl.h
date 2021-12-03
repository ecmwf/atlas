/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/functionspace.h"
#include "atlas/util/Object.h"

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
}  // namespace atlas

namespace atlas {
namespace redistribution {
namespace detail {

/// \brief  Abstract base class for redistributor implementation.
class RedistributionImpl : public util::Object {
public:
    RedistributionImpl() = default;

    /// \brief    Initialises the redistributor.
    ///
    /// \details  Performs MPI_Allgatherv to determine the (i, j, k) ranges
    ///           of each source and target function space on each PE.
    ///           The grids of source and target function space must match.
    ///
    /// \param[in]  source  Function space of source fields.
    /// \param[in]  target  Function space of target fields.
    void setup(const FunctionSpace& source, const FunctionSpace& target);

    /// \Setup class.
    virtual void do_setup() = 0;

    /// \brief Concrete type.
    virtual std::string type() const = 0;

    /// \brief  Maps source field to target field.
    virtual void execute(const Field& source, Field& target) const = 0;

    /// \brief  Maps source field set to target field set.
    virtual void execute(const FieldSet& source, FieldSet& target) const = 0;

    /// \brief  Get const reference to source function space.
    const FunctionSpace& source() const;

    /// \brief  Get const reference to target function space.
    const FunctionSpace& target() const;

private:
    FunctionSpace source_;
    FunctionSpace target_;
};

}  // namespace detail
}  // namespace redistribution
}  // namespace atlas
