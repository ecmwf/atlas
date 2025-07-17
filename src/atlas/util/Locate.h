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
#include <functional>

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/mdspan.h"

namespace atlas {

namespace detail {
struct DistributionConstView {
size_t size() const { return size_; }
int operator[](gidx_t i) const {
    return (*func_)(i);
}
DistributionConstView(std::function<int(gidx_t)>* func, std::size_t size) : func_(func), size_(size) {}
private:
std::function<int(gidx_t)>* func_{nullptr};
std::size_t size_{0};
};

#if ATLAS_HAVE_FORTRAN
constexpr size_t remote_index_base = 1;
#else
constexpr size_t remote_index_base = 0;
#endif
constexpr size_t global_index_base = 1;

template<class T, std::size_t Extents = dynamic_extent, class Layout = layout_right, class Accessor = default_accessor<T>>
using mdspan_1d = mdspan<T, extents<size_t,Extents>, Layout, Accessor>;

template<class T, std::size_t index_base , std::size_t Extents = dynamic_extent>
using mdspan_1d_index = mdspan_1d<T, Extents, layout_right, index_accessor<T,index_base>>;
}

namespace view {
using global_index = detail::mdspan_1d_index<gidx_t, detail::global_index_base>;
using remote_index = detail::mdspan_1d_index<idx_t,  detail::remote_index_base>;
using partition    = detail::mdspan_1d<int>;
}
namespace const_view {
using distribution = detail::DistributionConstView;
using global_index = detail::mdspan_1d_index<const gidx_t, detail::global_index_base>;
using remote_index = detail::mdspan_1d_index<const idx_t,  detail::remote_index_base>;
using partition    = detail::mdspan_1d<const int>;
using ghost        = detail::mdspan_1d<const int>;
}
}

namespace atlas::util {
// Given distribution and global_index, find the corresponding partition for each global index
// The global-index is typically 1-based
void locate_partition(
    // context
    const_view::distribution,
    // input
    const_view::global_index,
    // output
    view::partition);

// Given global_index, find the corresponding partition and remote_index
// This could be costly as it involves memory and communication
void locate_remote_index(
    // context
    std::string_view mpi_comm,
    const_view::global_index my_glb_idx, const_view::ghost my_ghost,
    // input
    const_view::global_index, const_view::partition,
    // output
    view::remote_index
);

// Given global_index, find the corresponding partition and remote_index
// This could be costly as it involves memory and communication
// local information: my_size, my_glb_idx, my_ghost
// global information: distribution
// requested indices to locate: size, global_index
// output of locate: partition, remote_index, remote_index_base
void locate(
    // context
    std::string_view mpi_comm,
    const_view::global_index my_glb_idx, const_view::ghost my_ghost,
    const_view::distribution,
    // input
    const_view::global_index,
    // output
    view::partition, view::remote_index);

// convenience functions
void locate(
    // context
    const FunctionSpace& fs,
    const_view::distribution,
    // input
    const_view::global_index,
    // output
    view::partition, view::remote_index);

void locate(
    // context
    const FunctionSpace& fs,
    // input
    const_view::global_index,
    // output
    view::partition, view::remote_index);

} // namespace atlas::util
