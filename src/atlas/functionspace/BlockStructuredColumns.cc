/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace/BlockStructuredColumns.h"

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------------

BlockStructuredColumns::BlockStructuredColumns(): FunctionSpace(), functionspace_(nullptr) {}

BlockStructuredColumns::BlockStructuredColumns(const FunctionSpace& functionspace):
    FunctionSpace(functionspace), 
    functionspace_(dynamic_cast<const detail::BlockStructuredColumns*>(get())) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const eckit::Configuration& config):
    FunctionSpace(new detail::BlockStructuredColumns(grid, config)),
    functionspace_(dynamic_cast<const detail::BlockStructuredColumns*>(get())) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Partitioner& partitioner,
                                     const eckit::Configuration& config):
    FunctionSpace(new detail::BlockStructuredColumns(grid, partitioner, config)),
    functionspace_(dynamic_cast<const detail::BlockStructuredColumns*>(get())) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Distribution& distribution,
                                     const eckit::Configuration& config):
    FunctionSpace(new detail::BlockStructuredColumns(grid, distribution, config)),
    functionspace_(dynamic_cast<const detail::BlockStructuredColumns*>(get())) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const Vertical& vertical, const eckit::Configuration& config):
    FunctionSpace(new detail::BlockStructuredColumns(grid, vertical, config)),
    functionspace_(dynamic_cast<const detail::BlockStructuredColumns*>(get())) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const Vertical& vertical, const grid::Partitioner& partitioner,
                                     const eckit::Configuration& config):
    FunctionSpace(new detail::BlockStructuredColumns(grid, vertical, partitioner, config)),
    functionspace_(dynamic_cast<const detail::BlockStructuredColumns*>(get())) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Distribution& distribution, const Vertical& vertical,
                                     const eckit::Configuration& config):
    FunctionSpace(new detail::BlockStructuredColumns(grid, distribution, vertical, config)),
    functionspace_(dynamic_cast<const detail::BlockStructuredColumns*>(get())) {}

std::string BlockStructuredColumns::checksum(const FieldSet& fieldset) const {
    return functionspace_->checksum(fieldset);
}

std::string BlockStructuredColumns::checksum(const Field& field) const {
    return functionspace_->checksum(field);
}

// ----------------------------------------------------------------------------

}  // namespace functionspace
}  // namespace atlas
