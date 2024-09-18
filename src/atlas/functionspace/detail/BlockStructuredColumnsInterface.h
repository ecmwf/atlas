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

#include "atlas/functionspace/BlockStructuredColumns.h"


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

const detail::BlockStructuredColumns* atlas__functionspace__BStructuredColumns__new__grid(const GridImpl* grid,
                                                                                    const eckit::Configuration* config);

const detail::BlockStructuredColumns* atlas__functionspace__BStructuredColumns__new__grid_dist(
    const GridImpl* grid, const grid::DistributionImpl* dist, const eckit::Configuration* config);

const detail::BlockStructuredColumns* atlas__functionspace__BStructuredColumns__new__grid_dist_vert(
    const GridImpl* grid, const grid::DistributionImpl* dist, const Vertical* vert, const eckit::Configuration* config);
const detail::BlockStructuredColumns* atlas__functionspace__BStructuredColumns__new__grid_dist_config(
    const GridImpl* grid, const grid::DistributionImpl* dist, const eckit::Configuration* config);

const detail::BlockStructuredColumns* atlas__functionspace__BStructuredColumns__new__grid_part(
    const GridImpl* grid, const PartitionerImpl* dist, const eckit::Configuration* config);

const detail::BlockStructuredColumns* atlas__functionspace__BStructuredColumns__new__grid_part_vert(
    const GridImpl* grid, const PartitionerImpl* dist, const Vertical* vert, const eckit::Configuration* config);

void atlas__functionspace__BStructuredColumns__delete(detail::BlockStructuredColumns* This);
field::FieldImpl* atlas__fs__BStructuredColumns__create_field(const detail::BlockStructuredColumns* This,
                                                             const eckit::Configuration* options);

void atlas__functionspace__BStructuredColumns__gather_field(const detail::BlockStructuredColumns* This,
                                                           const field::FieldImpl* local, field::FieldImpl* global);
void atlas__functionspace__BStructuredColumns__scatter_field(const detail::BlockStructuredColumns* This,
                                                            const field::FieldImpl* global, field::FieldImpl* local);
void atlas__functionspace__BStructuredColumns__gather_fieldset(const detail::BlockStructuredColumns* This,
                                                              const field::FieldSetImpl* local,
                                                              field::FieldSetImpl* global);
void atlas__functionspace__BStructuredColumns__scatter_fieldset(const detail::BlockStructuredColumns* This,
                                                               const field::FieldSetImpl* global,
                                                               field::FieldSetImpl* local);
void atlas__fs__BStructuredColumns__checksum_fieldset(const detail::BlockStructuredColumns* This,
                                                     const field::FieldSetImpl* fieldset, char*& checksum, idx_t& size,
                                                     int& allocated);
void atlas__fs__BStructuredColumns__checksum_field(const detail::BlockStructuredColumns* This, const field::FieldImpl* field,
                                                  char*& checksum, idx_t& size, int& allocated);
void atlas__fs__BStructuredColumns__index_host(const detail::BlockStructuredColumns* This, idx_t*& data, idx_t& i_min,
                                              idx_t& i_max, idx_t& j_min, idx_t& j_max);
idx_t atlas__fs__BStructuredColumns__size(const detail::BlockStructuredColumns* This);
idx_t atlas__fs__BStructuredColumns__sizeOwned(const detail::BlockStructuredColumns* This);
idx_t atlas__fs__BStructuredColumns__j_begin(const detail::BlockStructuredColumns* This);
idx_t atlas__fs__BStructuredColumns__j_end(const detail::BlockStructuredColumns* This);
idx_t atlas__fs__BStructuredColumns__i_begin(const detail::BlockStructuredColumns* This, idx_t j);
idx_t atlas__fs__BStructuredColumns__i_end(const detail::BlockStructuredColumns* This, idx_t j);
idx_t atlas__fs__BStructuredColumns__j_begin_halo(const detail::BlockStructuredColumns* This);
idx_t atlas__fs__BStructuredColumns__j_end_halo(const detail::BlockStructuredColumns* This);
idx_t atlas__fs__BStructuredColumns__i_begin_halo(const detail::BlockStructuredColumns* This, idx_t j);
idx_t atlas__fs__BStructuredColumns__i_end_halo(const detail::BlockStructuredColumns* This, idx_t j);
idx_t atlas__fs__BStructuredColumns__levels(const detail::BlockStructuredColumns* This);
idx_t atlas__fs__BStructuredColumns__block_begin(const detail::BlockStructuredColumns* This, idx_t jblk);
idx_t atlas__fs__BStructuredColumns__block_size(const detail::BlockStructuredColumns* This, idx_t jblk);
idx_t atlas__fs__BStructuredColumns__nproma(const detail::BlockStructuredColumns* This);
idx_t atlas__fs__BStructuredColumns__nblks(const detail::BlockStructuredColumns* This);

field::FieldImpl* atlas__fs__BStructuredColumns__xy(const detail::BlockStructuredColumns* This);
field::FieldImpl* atlas__fs__BStructuredColumns__z(const detail::BlockStructuredColumns* This);
field::FieldImpl* atlas__fs__BStructuredColumns__partition(const detail::BlockStructuredColumns* This);
field::FieldImpl* atlas__fs__BStructuredColumns__global_index(const detail::BlockStructuredColumns* This);
field::FieldImpl* atlas__fs__BStructuredColumns__index_i(const detail::BlockStructuredColumns* This);
field::FieldImpl* atlas__fs__BStructuredColumns__index_j(const detail::BlockStructuredColumns* This);

const GridImpl* atlas__fs__BStructuredColumns__grid(const detail::BlockStructuredColumns* This);
}

}  // namespace functionspace
}  // namespace atlas
