/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <vector>

#include "atlas/grid/StructuredGrid.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

StructuredGrid::grid_t* reduced_gaussian(const std::vector<long>& nx);
StructuredGrid::grid_t* reduced_gaussian(const std::vector<long>& nx, const Projection&);
StructuredGrid::grid_t* reduced_gaussian(const std::vector<long>& nx, const Domain&);

StructuredGrid::grid_t* reduced_gaussian(const std::vector<int>& nx);
StructuredGrid::grid_t* reduced_gaussian(const std::vector<int>& nx, const Projection&);
StructuredGrid::grid_t* reduced_gaussian(const std::vector<int>& nx, const Domain&);

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
