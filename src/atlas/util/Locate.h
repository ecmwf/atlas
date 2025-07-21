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

#include <string>
#include <string_view>
#include <vector>
#include <functional>

#include "atlas/mdspan.h"
#include "atlas/util/vector.h"

namespace atlas::parallel {

class Locator {
public:
    template <typename T>
    struct fspan {
        using value_type = std::remove_const_t<T>;
        size_t size() const { return size_; }
        value_type operator[](size_t i) const {
            return (*func_)(i);
        }
        fspan(const std::function<value_type(size_t)>* func, std::size_t size) :
            func_(func), size_(size) {
            }
    private:
        const std::function<value_type(size_t)>* func_{nullptr};
        std::size_t size_{0};
    };

    template<class T, std::size_t Extents = dynamic_extent, class Layout = layout_right, class Accessor = default_accessor<T>>
    using span = mdspan<T, extents<size_t,Extents>, Layout, Accessor>;

    static void locate_partition(
        // context
        fspan<const int> distribution,
        // input
        span<const gidx_t> global_index, const gidx_t global_index_base,
        // output
        span<int> partition);

    static void locate_partition(
        // context
        span<const int> distribution,
        // input
        span<const gidx_t> global_index, const gidx_t global_index_base,
        // output
        span<int> partition);

    // Given global_index, find the corresponding partition and remote_index
    // This could be costly as it involves memory and communication
    static void locate_remote_index(
        // context
        std::string_view mpi_comm,
        span<const gidx_t> my_glb_idx, const gidx_t my_global_index_base, span<const int> my_ghost,
        // input
        span<const gidx_t> global_index, const gidx_t global_index_base, span<const int> partition,
        // output
        span<idx_t> remote_index, const idx_t remote_index_base
    );

    // Given global_index, find the corresponding partition and remote_index
    // This could be costly as it involves memory and communication
    // local information: my_size, my_glb_idx, my_ghost
    // global information: distribution
    // requested indices to locate: size, global_index
    // output of locate: partition, remote_index, remote_index_base
    static void locate(
        // context
        std::string_view mpi_comm,
        span<const gidx_t> my_glb_idx, const gidx_t my_global_index_base, span<const int> my_ghost,
        fspan<const int> distribution,
        // input
        span<const gidx_t> global_index, const gidx_t global_index_base,
        // output
        span<int> partition, span<idx_t> remote_index, const idx_t remote_index_base);

    static void locate(
        // context
        std::string_view mpi_comm,
        span<const gidx_t> my_glb_idx, const gidx_t my_global_index_base, span<const int> my_ghost,
        span<const int> distribution,
        // input
        span<const gidx_t> global_index, const gidx_t global_index_base,
        // output
        span<int> partition, span<idx_t> remote_index, const idx_t remote_index_base);

    virtual void locate(
        // input
        span<const gidx_t> global_index, const gidx_t global_index_base,
        // output
        span<int> partition, span<idx_t> remote_index, const idx_t remote_index_base) const = 0;

    virtual void locate(
        // input
        const std::vector<gidx_t> global_index, const gidx_t global_index_base,
        // output
        std::vector<int>& partition, std::vector<idx_t>& remote_index, const idx_t remote_index_base) const {
            locate(
                span<const gidx_t>{global_index.data(), global_index.size()}, global_index_base,
                span<int>{partition.data(), partition.size()},
                span<idx_t>{remote_index.data(), remote_index.size()}, remote_index_base
            );
        }

};
} // namespace atlas::parallel

#include "atlas/functionspace/FunctionSpace.h"

namespace atlas::functionspace {
class Locator : public ::atlas::parallel::Locator {
public:
    Locator(FunctionSpace fs);

    using ::atlas::parallel::Locator::locate;
    void locate(
        // input
        span<const gidx_t> global_index, const gidx_t global_index_base,
        // output
        span<int> partition, span<idx_t> remote_index, const idx_t remote_index_base) const override;

private:
    FunctionSpace fs_;
    ::atlas::vector<int> distribution_array_;
    std::function<int(size_t)> distribution_function_;
    size_t distribution_size_;
};

} // namespace atlas::functionspace
