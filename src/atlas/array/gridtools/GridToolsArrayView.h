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

#include <cassert>
#include <cstddef>
#include <cstring>
#include <type_traits>

#include "atlas/array/Array.h"
#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/LocalView.h"
#include "atlas/array/gridtools/GridToolsMakeView.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

template <typename Value, int Rank>
class ArrayView {
    template <typename T>
    using is_non_const_value_type = typename std::is_same<T, typename std::remove_const<Value>::type>;

#define ENABLE_IF_NON_CONST                                                                             \
    template <bool EnableBool                                                                   = true, \
              typename std::enable_if<(!std::is_const<Value>::value && EnableBool), int>::type* = nullptr>

#define ENABLE_IF_CONST_WITH_NON_CONST(T)                                                                             \
    template <typename T, typename std::enable_if<(std::is_const<Value>::value && is_non_const_value_type<T>::value), \
                                                  int>::type* = nullptr>

public:
    // -- Type definitions
    using value_type                   = Value;
    using non_const_value_type         = typename std::remove_const<Value>::type;
    static constexpr bool is_const     = std::is_const<Value>::value;
    static constexpr bool is_non_const = !std::is_const<Value>::value;

    //    static constexpr Intent ACCESS{AccessMode};
    static constexpr int RANK{Rank};


    using data_view_t = gridtools::data_view_tt<Value, Rank, gridtools::get_access_mode<Value>()>;

private:
    using slicer_t       = typename helpers::ArraySlicer<ArrayView<Value, Rank>>;
    using const_slicer_t = typename helpers::ArraySlicer<const ArrayView<const Value, Rank>>;

    template <typename... Args>
    struct slice_t {
        using type = typename slicer_t::template Slice<Args...>::type;
    };

    template <typename... Args>
    struct const_slice_t {
        using type = typename const_slicer_t::template Slice<Args...>::type;
    };

public:
    ATLAS_HOST_DEVICE
    ArrayView(const ArrayView& other);
    ArrayView(const Array& array, bool device_view);

    ENABLE_IF_CONST_WITH_NON_CONST(value_type)
    ArrayView(const ArrayView<value_type, Rank>& other):
        gt_data_view_(other.is_device_view_ ? gridtools::make_gt_device_view<Value, Rank>(*other.array_)
                                            : gridtools::make_gt_host_view<Value, Rank>(*other.array_)),
        data_store_orig_(other.data_store_orig_),
        array_(other.array_),
        is_device_view_(other.is_device_view_) {
        std::memcpy(shape_, other.shape_, sizeof(ArrayShape::value_type) * Rank);
        std::memcpy(strides_, other.strides_, sizeof(ArrayStrides::value_type) * Rank);
        size_ = other.size_;
        // TODO: check compatibility
    }

    ENABLE_IF_CONST_WITH_NON_CONST(value_type)
    operator const ArrayView<value_type, Rank>&() const { return *(const ArrayView<value_type, Rank>*)(this); }

    value_type* data() { return gt_data_view_.data(); }
    value_type const* data() const { return gt_data_view_.data(); }

    template <typename... Ints, typename = typename std::enable_if<(sizeof...(Ints) == Rank), int>::type>
    ATLAS_HOST_DEVICE value_type& operator()(Ints... c) {
        assert(sizeof...(Ints) == Rank);
        using common_type = typename std::common_type<Ints...>::type;
        return gt_data_view_(static_cast<common_type>(c)...);
    }

    template <typename... Ints, typename = typename std::enable_if<(sizeof...(Ints) == Rank), int>::type>
    ATLAS_HOST_DEVICE value_type const& operator()(Ints... c) const {
        assert(sizeof...(Ints) == Rank);
        using common_type = typename std::common_type<Ints...>::type;
        return gt_data_view_(static_cast<common_type>(c)...);
    }

    template <typename Int, bool EnableBool = true>
    ATLAS_HOST_DEVICE typename std::enable_if<(Rank == 1 && EnableBool), const value_type&>::type operator[](
        Int idx) const {
        return gt_data_view_(idx);
    }

    template <typename Int, bool EnableBool = true>
    ATLAS_HOST_DEVICE typename std::enable_if<(Rank == 1 && EnableBool), value_type&>::type operator[](Int idx) {
        return gt_data_view_(idx);
    }

    template <unsigned int Dim>
    ATLAS_HOST_DEVICE idx_t shape() const {
        return gt_data_view_.template length<Dim>();
    }

    ATLAS_HOST_DEVICE
    data_view_t& data_view() { return gt_data_view_; }
    ATLAS_HOST_DEVICE
    data_view_t const& data_view() const { return gt_data_view_; }

    template <unsigned int Dim>
    ATLAS_HOST_DEVICE idx_t stride() const {
        return gt_data_view_.storage_info().template stride<Dim>();
    }

    static constexpr idx_t rank() { return Rank; }

    ATLAS_HOST_DEVICE size_t size() const { return size_; }
    bool valid() const;

    bool contiguous() const { return (size_ == size_t(shape_[0]) * size_t(strides_[0]) ? true : false); }

    void dump(std::ostream& os) const;

    ENABLE_IF_NON_CONST
    void assign(const value_type& value);

    ENABLE_IF_NON_CONST
    void assign(const std::initializer_list<value_type>& list);

    ENABLE_IF_NON_CONST
    void assign(const ArrayView<Value, Rank>& other);

    const idx_t* strides() const { return strides_; }

    const idx_t* shape() const { return shape_; }

    template <typename Int>
    ATLAS_HOST_DEVICE
    idx_t shape(Int idx) const {
        return shape_[idx];
    }

    template <typename Int>
    ATLAS_HOST_DEVICE
    idx_t stride(Int idx) const {
        return strides_[idx];
    }

    template <typename... Args>
    typename slice_t<Args...>::type slice(Args... args) {
        return slicer_t(*this).apply(args...);
    }

    template <typename... Args>
    typename const_slice_t<Args...>::type slice(Args... args) const {
        return const_slicer_t(*this).apply(args...);
    }

    bool isDeviceView() const { return is_device_view_; }

    // Befriend all template variations
    template <typename friendValue, int friendRank>
    friend class ArrayView;

private:
    data_view_t gt_data_view_;
    idx_t shape_[Rank];
    idx_t strides_[Rank];
    size_t size_;
    ArrayDataStore const* data_store_orig_;
    Array const* array_;
    bool is_device_view_{false};
#undef ENABLE_IF_NON_CONST
#undef ENABLE_IF_CONST_WITH_NON_CONST
};

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
