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
#include "atlas/util/mdspan.h"

namespace atlas::util {

template <typename Value, Value Base>
class IndexReference {
public:
    IndexReference(Value& idx): idx_(idx) {}

    void set(const Value& value) { idx_ = value + base_; }
    Value get() const { return idx_-base_; }
    void operator =(const Value& value) { set(value); }

    IndexReference& operator=(const IndexReference& other) {
        set(other.get());
        return *this;
    }
    IndexReference& operator--() {
        --(idx_);
        return *this;
    }
    IndexReference& operator++() {
        ++(idx_);
        return *this;
    }

    IndexReference& operator+=(Value v) {
        idx_ += v;
        return *this;
    }

    IndexReference& operator-=(Value v) {
        idx_ -= v;
        return *this;
    }

    // implicit conversion
    operator Value() const { return get(); }

private:
    Value& idx_;
    static constexpr Value base_ = Base;
};

template<class ElementType, size_t Base, typename Enable = void>
struct IndexAccessor;

template<class ElementType, size_t Base>
struct IndexAccessor<ElementType, Base, std::enable_if_t<!std::is_const_v<ElementType>>> {
    using element_type     = ElementType;
    using reference        = IndexReference<ElementType,Base>;
    using data_handle_type = ElementType*;
    using offset_policy    = IndexAccessor;

    constexpr IndexAccessor() = default;

    constexpr reference access(data_handle_type p, size_t i) const noexcept {
        return p[i];
    }

    constexpr data_handle_type offset(data_handle_type p, size_t i) const noexcept {
        return p + i;
    }

    constexpr size_t base() const noexcept { return base_; }

  private:
    static constexpr size_t base_{Base};                              // exposition only
};

template<class ElementType, size_t B>
struct IndexAccessor<ElementType, B, std::enable_if_t<std::is_const_v<ElementType>>>  {
    using element_type     = ElementType;
    using reference        = ElementType;
    using data_handle_type = ElementType*;
    using offset_policy    = IndexAccessor;

    constexpr IndexAccessor() = default;

    constexpr reference access(data_handle_type p, size_t i) const noexcept {
        return p[i] - base_;
    }

    constexpr data_handle_type offset(data_handle_type p, size_t i) const noexcept {
        return p + i;
    }

    constexpr size_t base() const noexcept { return base_; }

  private:
    static constexpr size_t base_{B};                              // exposition only
};

// Custom accessor policy that adds a base to the 0-based value
template<typename ElementType, typename Enable>
struct IndexAccessor<ElementType, 0, Enable> {
    using element_type     = ElementType;
    using reference        = ElementType&;
    using data_handle_type = ElementType*;
    using offset_policy    = IndexAccessor;

    constexpr IndexAccessor() = default;

    constexpr reference access(data_handle_type p, size_t i) const {
        return offset(p,i);
    }

    constexpr data_handle_type offset(data_handle_type p, size_t i) const noexcept {
        return p + i;
    }

    constexpr size_t base() const noexcept { return base_; }

  private:
    static constexpr size_t base_{0};                              // exposition only
};

template<class T, std::size_t Extents = dynamic_extent, class Accessor = default_accessor<T>>
using span = mdspan<T, extents<size_t,Extents>, layout_right, Accessor>;

template <size_t Base>
struct index_base {
    static constexpr size_t value = Base;
};

template<class T, class Base, std::size_t Extents = dynamic_extent>
using span_index = mdspan<T, extents<size_t,Extents>, layout_right, IndexAccessor<T,Base::value>>;

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

namespace view {
using global_index = span_index<gidx_t,index_base<1>>;
using remote_index = span_index<idx_t,index_base<1>>;
using partition    = span<int>;
}
namespace const_view {
using global_index = span_index<const gidx_t,index_base<1>>;
using remote_index = span_index<const idx_t,index_base<1>>;
using partition    = span<const int>;
using ghost        = span<const int>;
using distribution = DistributionConstView;
}


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
