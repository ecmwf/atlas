/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include "atlas/functionspace/StructuredColumns.h"


namespace atlas {
namespace field {
class FieldSetImpl;
}
namespace grid {
class DistributionImpl;
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
namespace grid {
namespace detail {
namespace partitioner {
class Partitioner;
}  // namespace partitioner
}  // namespace detail
}  // namespace grid
using PartitionerImpl = grid::detail::partitioner::Partitioner;
}  // namespace atlas


namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid(const GridImpl* grid,
                                                                                    const eckit::Configuration* config);

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_dist(
    const GridImpl* grid, const grid::DistributionImpl* dist, const eckit::Configuration* config);

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_dist_vert(
    const GridImpl* grid, const grid::DistributionImpl* dist, const Vertical* vert, const eckit::Configuration* config);
const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_dist_config(
    const GridImpl* grid, const grid::DistributionImpl* dist, const eckit::Configuration* config);

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_part(
    const GridImpl* grid, const PartitionerImpl* dist, const eckit::Configuration* config);

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_part_vert(
    const GridImpl* grid, const PartitionerImpl* dist, const Vertical* vert, const eckit::Configuration* config);

void atlas__functionspace__StructuredColumns__delete(detail::StructuredColumns* This);
field::FieldImpl* atlas__fs__StructuredColumns__create_field(const detail::StructuredColumns* This,
                                                             const eckit::Configuration* options);

void atlas__functionspace__StructuredColumns__gather_field(const detail::StructuredColumns* This,
                                                           const field::FieldImpl* local, field::FieldImpl* global);
void atlas__functionspace__StructuredColumns__scatter_field(const detail::StructuredColumns* This,
                                                            const field::FieldImpl* global, field::FieldImpl* local);
void atlas__functionspace__StructuredColumns__gather_fieldset(const detail::StructuredColumns* This,
                                                              const field::FieldSetImpl* local,
                                                              field::FieldSetImpl* global);
void atlas__functionspace__StructuredColumns__scatter_fieldset(const detail::StructuredColumns* This,
                                                               const field::FieldSetImpl* global,
                                                               field::FieldSetImpl* local);
void atlas__fs__StructuredColumns__checksum_fieldset(const detail::StructuredColumns* This,
                                                     const field::FieldSetImpl* fieldset, char*& checksum, idx_t& size,
                                                     int& allocated);
void atlas__fs__StructuredColumns__checksum_field(const detail::StructuredColumns* This, const field::FieldImpl* field,
                                                  char*& checksum, idx_t& size, int& allocated);
void atlas__fs__StructuredColumns__index_host(const detail::StructuredColumns* This, idx_t*& data, idx_t& i_min,
                                              idx_t& i_max, idx_t& j_min, idx_t& j_max);
idx_t atlas__fs__StructuredColumns__size(const detail::StructuredColumns* This);
idx_t atlas__fs__StructuredColumns__sizeOwned(const detail::StructuredColumns* This);
idx_t atlas__fs__StructuredColumns__j_begin(const detail::StructuredColumns* This);
idx_t atlas__fs__StructuredColumns__j_end(const detail::StructuredColumns* This);
idx_t atlas__fs__StructuredColumns__i_begin(const detail::StructuredColumns* This, idx_t j);
idx_t atlas__fs__StructuredColumns__i_end(const detail::StructuredColumns* This, idx_t j);
idx_t atlas__fs__StructuredColumns__j_begin_halo(const detail::StructuredColumns* This);
idx_t atlas__fs__StructuredColumns__j_end_halo(const detail::StructuredColumns* This);
idx_t atlas__fs__StructuredColumns__i_begin_halo(const detail::StructuredColumns* This, idx_t j);
idx_t atlas__fs__StructuredColumns__i_end_halo(const detail::StructuredColumns* This, idx_t j);
idx_t atlas__fs__StructuredColumns__levels(const detail::StructuredColumns* This);

field::FieldImpl* atlas__fs__StructuredColumns__xy(const detail::StructuredColumns* This);
field::FieldImpl* atlas__fs__StructuredColumns__z(const detail::StructuredColumns* This);
field::FieldImpl* atlas__fs__StructuredColumns__partition(const detail::StructuredColumns* This);
field::FieldImpl* atlas__fs__StructuredColumns__global_index(const detail::StructuredColumns* This);
field::FieldImpl* atlas__fs__StructuredColumns__index_i(const detail::StructuredColumns* This);
field::FieldImpl* atlas__fs__StructuredColumns__index_j(const detail::StructuredColumns* This);

const GridImpl* atlas__fs__StructuredColumns__grid(const detail::StructuredColumns* This);
}

}  // namespace functionspace
}  // namespace atlas
