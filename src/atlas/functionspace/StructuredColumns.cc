/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace/StructuredColumns.h"

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------------

StructuredColumns::StructuredColumns(): FunctionSpace(), functionspace_(nullptr) {}

StructuredColumns::StructuredColumns(const FunctionSpace& functionspace):
    FunctionSpace(functionspace), functionspace_(dynamic_cast<const detail::StructuredColumns*>(get())) {}

StructuredColumns::StructuredColumns(const Grid& grid, const eckit::Configuration& config):
    FunctionSpace(new detail::StructuredColumns(grid, config)),
    functionspace_(dynamic_cast<const detail::StructuredColumns*>(get())) {}

StructuredColumns::StructuredColumns(const Grid& grid, const grid::Partitioner& partitioner,
                                     const eckit::Configuration& config):
    FunctionSpace(new detail::StructuredColumns(grid, partitioner, config)),
    functionspace_(dynamic_cast<const detail::StructuredColumns*>(get())) {}

StructuredColumns::StructuredColumns(const Grid& grid, const grid::Distribution& distribution,
                                     const eckit::Configuration& config):
    FunctionSpace(new detail::StructuredColumns(grid, distribution, config)),
    functionspace_(dynamic_cast<const detail::StructuredColumns*>(get())) {}

StructuredColumns::StructuredColumns(const Grid& grid, const Vertical& vertical, const eckit::Configuration& config):
    FunctionSpace(new detail::StructuredColumns(grid, vertical, config)),
    functionspace_(dynamic_cast<const detail::StructuredColumns*>(get())) {}

StructuredColumns::StructuredColumns(const Grid& grid, const Vertical& vertical, const grid::Partitioner& partitioner,
                                     const eckit::Configuration& config):
    FunctionSpace(new detail::StructuredColumns(grid, vertical, partitioner, config)),
    functionspace_(dynamic_cast<const detail::StructuredColumns*>(get())) {}

StructuredColumns::StructuredColumns(const Grid& grid, const grid::Distribution& distribution, const Vertical& vertical,
                                     const eckit::Configuration& config):
    FunctionSpace(new detail::StructuredColumns(grid, distribution, vertical, config)),
    functionspace_(dynamic_cast<const detail::StructuredColumns*>(get())) {}

std::string StructuredColumns::checksum(const FieldSet& fieldset) const {
    return functionspace_->checksum(fieldset);
}

std::string StructuredColumns::checksum(const Field& field) const {
    return functionspace_->checksum(field);
}

// ----------------------------------------------------------------------------

}  // namespace functionspace
}  // namespace atlas
