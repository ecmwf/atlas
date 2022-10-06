/*
 * (C) Copyright 200 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include "atlas/functionspace/PointCloud.h"


namespace atlas {
namespace field {
class FieldSetImpl;
}
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
class Grid;
}  // namespace grid
}  // namespace detail
}  // namespace grid
using GridImpl = grid::detail::grid::Grid;
}  // namespace atlas

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {

const detail::PointCloud* atlas__functionspace__PointCloud__new__lonlat(const field::FieldImpl* lonlat);

const detail::PointCloud* atlas__functionspace__PointCloud__new__lonlat_ghost(const field::FieldImpl* lonlat,
                                                                              const field::FieldImpl* ghost);

const detail::PointCloud* atlas__functionspace__PointCloud__new__fieldset(const field::FieldSetImpl* fset);

const detail::PointCloud* atlas__functionspace__PointCloud__new__grid(const GridImpl* grid);

void atlas__functionspace__PointCloud__delete(detail::PointCloud* This);
field::FieldImpl* atlas__fs__PointCloud__create_field(const detail::PointCloud* This,
                                                      const eckit::Configuration* options);

idx_t atlas__fs__PointCloud__size(const detail::PointCloud* This);

const field::FieldImpl* atlas__fs__PointCloud__lonlat(const detail::PointCloud* This);
}

}  // namespace functionspace
}  // namespace atlas
