/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/runtime/Log.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/Range.h"

namespace atlas {
namespace array {

template< typename Value, int Rank, Intent AccessMode >
class LocalView;


template<typename Value>
struct Reference {
  Value& value;
  operator Value&() { return value; }
  template < typename T > void operator=(const T a) { value = a; }
  template < typename T > Reference<Value>& operator+(const T a) { value+=a; return *this; }
  template < typename T > Reference<Value>& operator-(const T a) { value-=a; return *this; }
  Reference<Value>& operator--() { --value; return *this; }
  Reference<Value>& operator++() { ++value; return *this; }
};


template<typename Value, int Rank, Intent Access>
struct get_slice_type {
  using type =
    typename std::conditional<(Rank==0),
       typename std::conditional<(Access==Intent::ReadOnly),
           Reference<Value const>,
           Reference<Value>
       >::type,
       LocalView<Value,Rank,Access>
    >::type;
};

//------------------------------------------------------------------------------

namespace helpers {

template< int Dim >
struct deduce_slice_rank;

template< int Rank, typename ...Args >
struct SliceRank_impl;

template<>
struct deduce_slice_rank<1> {

  template<typename Last>
  static constexpr int apply() {
    return std::is_base_of<RangeBase,Last>::value;
  }

};
template<typename ...Args>
struct SliceRank_impl<1,Args...> {
  static constexpr int value{ deduce_slice_rank<1>::apply<Args...>() };
};


#define ATLAS_ARRAY_SLICER_EXPLICIT_TEMPLATE_SPECIALISATION(RANK) \
template<> \
struct deduce_slice_rank<RANK> { \
\
  template<typename First, typename ...Args> \
  static constexpr int apply() { \
    return std::is_base_of<RangeBase,First>::value + deduce_slice_rank<RANK-1>::apply<Args...>(); \
  } \
 \
}; \
template<typename ...Args> \
struct SliceRank_impl<RANK,Args...> \
{ \
  static constexpr int value{ deduce_slice_rank<RANK>::apply<Args...>() }; \
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



template< typename View, bool constness = false >
struct get_access {
  static constexpr Intent value{ View::ACCESS };
};

template< typename View >
struct get_access<View,true> {
  static constexpr Intent value{ Intent::ReadOnly };
};


//template <typename Value, int Rank, Intent AccessMode>
template <typename View>
class ArraySlicer {
public:

    ArraySlicer(View& view) :
      view_(view) {
    }

    template <typename ...Args>
    struct SliceRank {
      static constexpr int value{ SliceRank_impl<View::RANK,Args...>::value };
    };

    template <typename ...Args>
    struct Slice {
      using type = typename get_slice_type<
          typename View::value_type,
          SliceRank<Args...>::value,
          get_access<View, std::is_const<View>::value>::value >::type;
    };

    template <typename ...Args>
    typename Slice<Args...>::type apply(const Args ...args) const {
        using slicer_t = Slicer< typename Slice<Args...>::type, (SliceRank<Args...>::value==0) >;
        return slicer_t::apply( view_, args... );
    }

private:

    template <typename ...Args>
    struct array {
      using type = typename std::array<size_t, SliceRank<Args...>::value >;
    };

    template <typename ReturnType, bool ToScalar = false>
    struct Slicer {
        template< typename ...Args >
        static ReturnType apply(View& view, const Args...args) {
            return ReturnType(
                  view.data()+offset(view,args...),
                  shape(view,args...).data(),
                  strides(view,args...).data()
            );
        }
    };

    template <typename ReturnType>
    struct Slicer<ReturnType,true> {
      template< typename ...Args >
      static ReturnType apply(View& view, const Args...args) {
          return ReturnType{*(view.data()+offset(view,args...))};
      }
    };

    template< typename Int >
    static int start( Int idx ) {
      return idx;
    }

    static int start( Range range ) {
      return range.start();
    }

    static int start( RangeAll range ) {
      return range.start();
    }

    static int start( RangeTo range ) {
      return range.start();
    }

    static int start( RangeFrom range ) {
      return range.start();
    }

    template < int Dim, typename Int, typename... Ints >
    static constexpr int offset_part(View& view, const Int idx, const Ints... next_idx) {
        return start(idx)*view.stride(Dim) + offset_part<Dim+1>( view, next_idx... );
    }

    template < int Dim, typename Int >
    static constexpr int offset_part(View& view, const Int last_idx) {
        return start(last_idx)*view.stride(Dim);
    }

    template < typename... Args >
    static constexpr int offset(View& view, const Args... args) {
        return offset_part<0>(view,args...);
    }


    template< int Dim, typename Shape, typename Int >
    static void update_shape( View&, Shape&, int& /*i*/, const Int& /*index*/ ){
      // do nothing
    }
    template< int Dim, typename Shape >
    static void update_shape( View&, Shape& shape, int& i, const Range range ){
       shape[i++] = range.end()-range.start();
    }
    template< int Dim, typename Shape >
    static void update_shape( View& view, Shape& shape, int& i, const RangeAll range ){
       shape[i++] = range.end<Dim>(view)-range.start();
    }
    template< int Dim, typename Shape >
    static void update_shape( View& view, Shape& shape, int& i, const RangeFrom range ){
       shape[i++] = range.end<Dim>(view)-range.start();
    }
    template< int Dim, typename Shape >
    static void update_shape( View&, Shape& shape, int& i, const RangeTo range ){
       shape[i++] = range.end()-range.start();
    }

    template< int Dim, typename Shape, typename Int, typename... Ints >
    static void shape_part( View& view, Shape& shape, int& i, const Int idx, const Ints... next_idx ) {
      update_shape<Dim>(view,shape,i,idx);
      shape_part<Dim+1>(view,shape,i,next_idx...);
    }

    template< int Dim, typename Shape, typename Int >
    static void shape_part( View& view, Shape& shape, int& i, const Int idx) {
      update_shape<Dim>(view,shape,i,idx);
    }

    template< typename... Args >
    static typename array<Args...>::type shape( View& view, const Args... args ) {
      typename array<Args...>::type result;
      int i(0);
      shape_part<0>(view,result,i,args...);
      return result;
    }


    template<int Dim, typename Strides, typename Int >
    static void update_strides( View&, Strides&, int& /*i*/, const Int& /*idx*/) {
      // do nothing
     }
    template<int Dim,typename Strides >
    static void update_strides( View& view, Strides& strides, int& i, const Range& /*range*/ ) {
       strides[i++] = view.stride(Dim);
    }
    template<int Dim,typename Strides >
    static void update_strides( View& view, Strides& strides, int& i, const RangeFrom& /*range*/ ) {
       strides[i++] = view.stride(Dim);
    }
    template<int Dim,typename Strides >
    static void update_strides( View& view, Strides& strides, int& i, const RangeTo& /*range*/ ) {
       strides[i++] = view.stride(Dim);
    }
    template<int Dim,typename Strides >
    static void update_strides( View& view, Strides& strides, int& i, const RangeAll& /*range*/ ) {
       strides[i++] = view.stride(Dim);
    }

    template< int Dim, typename Strides, typename Int, typename... Ints >
    static void strides_part( View& view, Strides& strides, int& i, const Int idx, const Ints... next_idx ) {
      update_strides<Dim>(view,strides,i,idx);
      strides_part<Dim+1>(view,strides,i,next_idx...);
    }

    template< int Dim, typename Strides, typename Int >
    static void strides_part( View& view, Strides& strides, int& i, const Int idx) {
      update_strides<Dim>(view,strides,i,idx);
    }

    template< typename... Args >
    static typename array<Args...>::type strides(View& view, const Args... args ) {
      typename array<Args...>::type result;
      int i(0);
      strides_part<0>(view,result,i,args...);
      return result;
    }

private:
    View& view_;
};

//------------------------------------------------------------------------------

} // namespace helpers
} // namespace array
} // namespace atlas
