/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas::util {


// Given distribution and global_index, find the corresponding partition for each global index
// The global-index is typically 1-based
void locate_partition(
    const grid::Distribution& distribution,
    std::size_t size, const gidx_t global_index[], gidx_t global_index_base,
    // output
    int partition[]);


// Given global_index, find the corresponding partition and remote_index
// This could be costly as it involves memory and communication
void locate_remote_index(
    const mpi::Comm& comm,
    std::size_t my_size, const gidx_t my_glb_idx[], const int my_ghost[],
    std::size_t size, const gidx_t global_index[], const int partition[],
    // output
    idx_t remote_index[], idx_t remote_index_base);

// Given global_index, find the corresponding partition and remote_index
// This could be costly as it involves memory and communication
// local information: my_size, my_glb_idx, my_ghost
// global information: distribution
// requested indices to locate: size, global_index
// output of locate: partition, remote_index, remote_index_base
void locate(const mpi::Comm& comm,
    std::size_t my_size, const gidx_t my_glb_idx[], int my_ghost[],
    const grid::Distribution& distribution,
    std::size_t size, const gidx_t global_index[], gidx_t global_index_base,
    // output
    int partition[], idx_t remote_index[], idx_t remote_index_base);


void locate(
    const FunctionSpace& fs, const grid::Distribution& distribution, std::size_t size,
    const gidx_t global_index[], gidx_t global_index_base,
    // output
    int partition[], idx_t remote_index[], idx_t remote_index_base);


void locate(const FunctionSpace& fs, std::size_t size, const gidx_t global_index[], const gidx_t global_index_base,
    // output
    int partition[], idx_t remote_index[], idx_t remote_index_base);

void locate(const FunctionSpace& fs, const std::vector<gidx_t>& global_index, gidx_t global_index_base,
    // output
    std::vector<int>& partition, std::vector<idx_t>& remote_index, idx_t remote_index_base);

void locate(const FunctionSpace& fs, const grid::Distribution& distribution, const std::vector<gidx_t>& global_index, gidx_t global_index_base,
    // output
    std::vector<int>& partition, std::vector<idx_t>& remote_index, idx_t remote_index_base);

} // namespace atlas::util
