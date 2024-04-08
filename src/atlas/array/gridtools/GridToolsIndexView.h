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
#include <type_traits>

#include "atlas/array/Array.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

namespace detail {

// FortranIndex:
// Helper class that does +1 and -1 operations on stored values
template <typename Value>
class FortranIndex {
public:
    enum
    {
        BASE = 1
    };

public:
    ATLAS_HOST_DEVICE
    FortranIndex(Value* idx): idx_(idx) {}
    ATLAS_HOST_DEVICE
    void set(const Value& value) { *(idx_) = value + BASE; }
    ATLAS_HOST_DEVICE
    Value get() const { return *(idx_)-BASE; }
    ATLAS_HOST_DEVICE
    void operator=(const Value& value) { set(value); }
    ATLAS_HOST_DEVICE
    FortranIndex<Value>& operator=(const FortranIndex<Value>& other) {
        set(other.get());
        return *this;
    }
    ATLAS_HOST_DEVICE
    FortranIndex<Value>& operator--() {
        --(*(idx_));
        return *this;
    }
    ATLAS_HOST_DEVICE
    FortranIndex<Value>& operator++() {
        ++(*(idx_));
        return *this;
    }

    ATLAS_HOST_DEVICE
    FortranIndex<Value>& operator+=(Value v) {
        *(idx_) += v;
        return *this;
    }

    ATLAS_HOST_DEVICE
    FortranIndex<Value>& operator-=(Value v) {
        *(idx_) -= v;
        return *this;
    }

    // implicit conversion
    ATLAS_HOST_DEVICE
    operator Value() const { return get(); }

private:
    Value* idx_;
};

}  // namespace detail

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
class IndexView {
public:
// -- Type definitions
#if ATLAS_HAVE_FORTRAN
    typedef detail::FortranIndex<Value> Index;
#define INDEX_REF Index
#define FROM_FORTRAN -1
#define TO_FORTRAN +1
#else
    typedef Value& Index;
#define INDEX_REF *
#define FROM_FORTRAN
#define TO_FORTRAN
#endif
    using data_view_t =
        gridtools::data_view_tt<typename std::remove_const<Value>::type, Rank, gridtools::get_access_mode<Value>()>;

    using value_type                   = Value;

public:
    IndexView(data_view_t);
    IndexView(const Array&, bool device_view);

    template <typename... Coords, typename = typename std::enable_if<(sizeof...(Coords) == Rank), int>::type>
    ATLAS_HOST_DEVICE
    Index operator()(Coords... c) {
        assert(sizeof...(Coords) == Rank);
        return INDEX_REF(&gt_data_view_(c...));
    }

    template <typename... Coords, typename = typename std::enable_if<(sizeof...(Coords) == Rank), int>::type>
    ATLAS_HOST_DEVICE
    Value const operator()(Coords... c) const {
        return gt_data_view_(c...) FROM_FORTRAN;
    }

    ATLAS_HOST_DEVICE
    idx_t size() const { return size_; }

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


    void dump(std::ostream& os) const;

private:
    data_view_t gt_data_view_;
    idx_t size_;
    idx_t shape_[Rank];
    idx_t strides_[Rank];
    bool is_device_view_;
    ArrayDataStore const* data_store_orig_;
    Array const* array_;

#undef INDEX_REF
#undef TO_FORTRAN
#undef FROM_FORTRAN
};


#if ATLAS_HAVE_FORTRAN
#define INDEX_REF Index
#define FROM_FORTRAN -1
#define TO_FORTRAN +1
#else
#define INDEX_REF *
#define FROM_FORTRAN
#define TO_FORTRAN
#endif


template <typename Value, int Rank>
class LocalIndexView {
public:
    using value_type = typename remove_const<Value>::type;

#if ATLAS_HAVE_FORTRAN
    typedef detail::FortranIndex<Value> Index;
#else
    typedef Value& Index;
#endif

public:
    LocalIndexView(Value* data, const idx_t shape[Rank]);

    LocalIndexView(Value* data, const idx_t shape[Rank], const idx_t strides[Rank]);

    // -- Access methods

    template <typename... Ints, typename = typename std::enable_if<(sizeof...(Ints) == Rank), int>::type>
    Index operator()(Ints... idx) {
        check_bounds(idx...);
        return INDEX_REF(&data_[index(idx...)]);
    }

    template <typename... Ints, typename = typename std::enable_if<(sizeof...(Ints) == Rank), int>::type>
    const value_type operator()(Ints... idx) const {
        return data_[index(idx...)] FROM_FORTRAN;
    }

private:
    // -- Private methods

    template <int Dim, typename Int, typename... Ints>
    constexpr idx_t index_part(Int idx, Ints... next_idx) const {
        return idx * strides_[Dim] + index_part<Dim + 1>(next_idx...);
    }

    template <int Dim, typename Int>
    constexpr idx_t index_part(Int last_idx) const {
        return last_idx * strides_[Dim];
    }

    template <typename... Ints>
    constexpr idx_t index(Ints... idx) const {
        return index_part<0>(idx...);
    }

#if ATLAS_INDEXVIEW_BOUNDS_CHECKING
    template <typename... Ints>
    void check_bounds(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
        return check_bounds_part<0>(idx...);
    }
#else
    template <typename... Ints>
    void check_bounds(Ints...) const {}
#endif

    template <typename... Ints>
    void check_bounds_force(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
        return check_bounds_part<0>(idx...);
    }

    template <int Dim, typename Int, typename... Ints>
    void check_bounds_part(Int idx, Ints... next_idx) const {
        if (idx_t(idx) >= shape_[Dim]) {
            throw_OutOfRange("IndexView", array_dim<Dim>(), idx, shape_[Dim]);
        }
        check_bounds_part<Dim + 1>(next_idx...);
    }

    template <int Dim, typename Int>
    void check_bounds_part(Int last_idx) const {
        if (idx_t(last_idx) >= shape_[Dim]) {
            throw_OutOfRange("IndexView", array_dim<Dim>(), last_idx, shape_[Dim]);
        }
    }

    idx_t size() const { return shape_[0]; }

    void dump(std::ostream& os) const;

private:
    Value* data_;
    idx_t strides_[Rank];
    idx_t shape_[Rank];
};

#undef INDEX_REF
#undef FROM_FORTRAN
#undef TO_FORTRAN

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
