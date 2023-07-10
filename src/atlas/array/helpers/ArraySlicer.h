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

#include <array>

#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/Range.h"
#include "atlas/library/config.h"

namespace atlas {
namespace array {

template <typename Value, int Rank>
class LocalView;

template <typename Value>
struct Reference {
    Value& value;
    operator Value&() { return value; }
    operator const Value&() const { return value; }
    template <typename T>
    void operator=(const T a) {
        value = a;
    }
    template <typename T>
    Reference<Value>& operator+(const T a) {
        value += a;
        return *this;
    }
    template <typename T>
    Reference<Value>& operator-(const T a) {
        value -= a;
        return *this;
    }
    Reference<Value>& operator--() {
        --value;
        return *this;
    }
    Reference<Value>& operator++() {
        ++value;
        return *this;
    }
    template <typename T>
    bool operator==(const T a) {
        return value == a;
    }
    template <typename T>
    bool operator!=(const T a) {
        return value != a;
    }
    friend std::ostream& operator<<(std::ostream& out, const Reference& ref) {
        out << ref.value;
        return out;
    }
    constexpr int shape(idx_t) { return 0; }
    constexpr int stride(idx_t) { return 0; }
    constexpr int rank() { return 0; }
};

template <typename Value, int Rank>
struct get_slice_type {
    using type = typename std::conditional<(Rank == 0), Reference<Value>, LocalView<Value, Rank>>::type;
};

//------------------------------------------------------------------------------

namespace helpers {

template <int Dim>
struct deduce_slice_rank;

template <int Rank, typename... Args>
struct SliceRank_impl;

template <>
struct deduce_slice_rank<1> {
    template <typename Last>
    static constexpr int apply() {
        return std::is_base_of<RangeBase, Last>::value;
    }
};
template <typename... Args>
struct SliceRank_impl<1, Args...> {
    static constexpr int value{deduce_slice_rank<1>::apply<Args...>()};
};

#define ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(RANK)                                            \
    template <>                                                                                              \
    struct deduce_slice_rank<RANK> {                                                                         \
        template <typename First, typename... Args>                                                          \
        static constexpr int apply() {                                                                       \
            return std::is_base_of<RangeBase, First>::value + deduce_slice_rank<RANK - 1>::apply<Args...>(); \
        }                                                                                                    \
    };                                                                                                       \
    template <typename... Args>                                                                              \
    struct SliceRank_impl<RANK, Args...> {                                                                   \
        static constexpr int value{deduce_slice_rank<RANK>::apply<Args...>()};                               \
    }

ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(2);
ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(3);
ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(4);
ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(5);
ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(6);
ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(7);
ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(8);
ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(9);
#undef ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION

template <typename... Args>
struct SliceRank {
    static constexpr int value{SliceRank_impl<sizeof...(Args), Args...>::value};
};

// template <typename Value, int Rank, Intent AccessMode>
template <typename View>
class ArraySlicer {
public:
    ArraySlicer(View& view): view_(view) {}

    template <typename... Args>
    struct Slice {
        using type = typename get_slice_type<typename View::value_type, SliceRank<Args...>::value>::type;
    };

    template <typename... Args>
    typename Slice<Args...>::type apply(const Args... args) const {
        using slicer_t = Slicer<typename Slice<Args...>::type, (SliceRank<Args...>::value == 0)>;
        static_assert(
            View::RANK <= sizeof...(Args),
            "Not enough arguments passed to slice() function. Pehaps you forgot to add a array::Range::all()");
        return slicer_t::apply(view_, args...);
    }

private:
    template <typename... Args>
    struct array {
        using type = typename std::array<idx_t, SliceRank<Args...>::value>;
    };

    template <typename ReturnType, bool ToScalar = false>
    struct Slicer {
        template <typename... Args>
        static ReturnType apply(View& view, const Args... args) {
            return ReturnType(view.data() + offset(view, args...), shape(view, args...).data(),
                              strides(view, args...).data());
        }
    };

    template <typename ReturnType>
    struct Slicer<ReturnType, true> {
        template <typename... Args>
        static ReturnType apply(View& view, const Args... args) {
            return ReturnType{*(view.data() + offset(view, args...))};
        }
    };

    template <typename Int>
    static idx_t offset_part(View& view, int& i_view, Int idx) {
        return idx * view.stride(i_view++);
    }

    static idx_t offset_part(View& view, int& i_view, Range range) { return range.start() * view.stride(i_view++); }

    static idx_t offset_part(View& view, int& i_view, RangeAll range) { return range.start() * view.stride(i_view++); }

    static idx_t offset_part(View& view, int& i_view, RangeTo range) { return range.start() * view.stride(i_view++); }

    static idx_t offset_part(View& view, int& i_view, RangeFrom range) { return range.start() * view.stride(i_view++); }

    static idx_t offset_part(View&, int& /*i_view*/, RangeDummy) { return 0; }

    template <int Dim, typename Int, typename... Ints>
    static idx_t offset_remaining(View& view, int& i_view, const Int idx, const Ints... next_idx) {
        const idx_t p = offset_part(view, i_view, idx);
        return p + offset_remaining<Dim + 1>(view, i_view, next_idx...);
    }

    template <int Dim, typename Int>
    static idx_t offset_remaining(View& view, int& i_view, const Int last_idx) {
        return offset_part(view, i_view, last_idx);
    }

    template <typename... Args>
    static idx_t offset(View& view, const Args... args) {
        int i_view(0);
        return offset_remaining<0>(view, i_view, args...);
    }

    template <int Dim, typename Shape, typename Int>
    static void update_shape(View&, Shape&, int& i_view, int& /*i_slice*/, const Int& /*index*/) {
        // do nothing
        ++i_view;
    }
    template <int Dim, typename Shape>
    static void update_shape(View&, Shape& shape, int& i_view, int& i_slice, const Range range) {
        shape[i_slice] = range.end() - range.start();
        ++i_slice;
        ++i_view;
    }
    template <int Dim, typename Shape>
    static void update_shape(View& view, Shape& shape, int& i_view, int& i_slice, const RangeAll range) {
        shape[i_slice] = range.end(view, i_view) - range.start();
        ++i_slice;
        ++i_view;
    }
    template <int Dim, typename Shape>
    static void update_shape(View& view, Shape& shape, int& i_view, int& i_slice, const RangeFrom range) {
        shape[i_slice] = range.end(view, i_view) - range.start();
        ++i_slice;
        ++i_view;
    }
    template <int Dim, typename Shape>
    static void update_shape(View&, Shape& shape, int& i_view, int& i_slice, const RangeTo range) {
        shape[i_slice] = range.end() - range.start();
        ++i_slice;
        ++i_view;
    }

    template <int Dim, typename Shape>
    static void update_shape(View&, Shape& shape, int& /*i_view*/, int& i_slice, const RangeDummy) {
        shape[i_slice] = 1;
        ++i_slice;
        // no update of i_view for dummy-dimension
    }

    template <int Dim, typename Shape, typename Int, typename... Ints>
    static void shape_part(View& view, Shape& shape, int& i_view, int& i_slice, const Int idx, const Ints... next_idx) {
        update_shape<Dim>(view, shape, i_view, i_slice, idx);
        shape_part<Dim + 1>(view, shape, i_view, i_slice, next_idx...);
    }

    template <int Dim, typename Shape, typename Int>
    static void shape_part(View& view, Shape& shape, int& i_view, int& i_slice, const Int idx) {
        update_shape<Dim>(view, shape, i_view, i_slice, idx);
    }

    template <typename... Args>
    static typename array<Args...>::type shape(View& view, const Args... args) {
        typename array<Args...>::type result;
        int i_slice(0);
        int i_view(0);
        shape_part<0>(view, result, i_view, i_slice, args...);
        return result;
    }

    template <int Dim, typename Strides, typename Int>
    static void update_strides(View&, Strides&, int& i_view, int& /*i_slice*/, const Int& /*idx*/) {
        // do nothing
        ++i_view;
    }
    template <int Dim, typename Strides>
    static void update_strides(View& view, Strides& strides, int& i_view, int& i_slice, const Range& /*range*/) {
        strides[i_slice] = view.stride(i_view);
        ++i_slice;
        ++i_view;
    }
    template <int Dim, typename Strides>
    static void update_strides(View& view, Strides& strides, int& i_view, int& i_slice, const RangeFrom& /*range*/) {
        strides[i_slice] = view.stride(i_view);
        ++i_slice;
        ++i_view;
    }
    template <int Dim, typename Strides>
    static void update_strides(View& view, Strides& strides, int& i_view, int& i_slice, const RangeTo& /*range*/) {
        strides[i_slice] = view.stride(i_view);
        ++i_slice;
        ++i_view;
    }
    template <int Dim, typename Strides>
    static void update_strides(View& view, Strides& strides, int& i_view, int& i_slice, const RangeAll& /*range*/) {
        strides[i_slice] = view.stride(i_view);
        ++i_slice;
        ++i_view;
    }
    template <int Dim, typename Strides>
    static void update_strides(View& /*view*/, Strides& strides, int& /*i_view*/, int& i_slice,
                               const RangeDummy& /*range*/) {
        strides[i_slice] = 0;
        ++i_slice;
    }

    template <int Dim, typename Strides, typename Int, typename... Ints>
    static void strides_part(View& view, Strides& strides, int& i_view, int& i_slice, const Int idx,
                             const Ints... next_idx) {
        update_strides<Dim>(view, strides, i_view, i_slice, idx);
        strides_part<Dim + 1>(view, strides, i_view, i_slice, next_idx...);
    }

    template <int Dim, typename Strides, typename Int>
    static void strides_part(View& view, Strides& strides, int& i_view, int& i_slice, const Int idx) {
        update_strides<Dim>(view, strides, i_view, i_slice, idx);
    }

    template <typename... Args>
    static typename array<Args...>::type strides(View& view, const Args... args) {
        typename array<Args...>::type result;
        int i_slice(0);
        int i_view(0);
        strides_part<0>(view, result, i_view, i_slice, args...);
        return result;
    }

private:
    View& view_;
};

//------------------------------------------------------------------------------

}  // namespace helpers
}  // namespace array
}  // namespace atlas
