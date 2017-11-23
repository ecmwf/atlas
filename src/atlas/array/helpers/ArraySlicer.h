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

#include <iostream>

#include "atlas/array.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace helpers {


class Range {
public:
  int start;
  int end;
};

class Index {
public:
  int index;
};

template< int Rank, typename ...Args >
struct SliceRank;




template< int Dim >
struct deduce_slice_rank;


template<>
struct deduce_slice_rank<1> {

  template<typename Last>
  static constexpr int apply() {
    return std::is_same<Last,Range>::value;
  }

};

#define EXPLICIT_TEMPLATE_SPECIALISATION(RANK) \
template<> \
struct deduce_slice_rank<RANK> { \
\
  template<typename First, typename ...Args> \
  static constexpr int apply() { \
    return std::is_same<First,Range>::value + deduce_slice_rank<RANK-1>::apply<Args...>(); \
  } \
 \
}; \
template<typename ...Args> \
struct SliceRank<RANK,Args...> \
{ \
  enum { value = deduce_slice_rank<RANK>::apply<Args...>() }; \
}

EXPLICIT_TEMPLATE_SPECIALISATION(2);
EXPLICIT_TEMPLATE_SPECIALISATION(3);
EXPLICIT_TEMPLATE_SPECIALISATION(4);
EXPLICIT_TEMPLATE_SPECIALISATION(5);
EXPLICIT_TEMPLATE_SPECIALISATION(6);
EXPLICIT_TEMPLATE_SPECIALISATION(7);
EXPLICIT_TEMPLATE_SPECIALISATION(8);
EXPLICIT_TEMPLATE_SPECIALISATION(9);
#undef EXPLICIT_TEMPLATE_SPECIALISATION

template<typename ...Args>
struct SliceRank<1,Args...>
{
  enum { value = deduce_slice_rank<1>::apply<Args...>() };
};




template <typename Value, int Rank, Intent AccessMode>
class ArraySlicer {
public:
    using value_type = Value;
    template<int RANK, Intent Access = AccessMode>
    struct Slice {
      using type = typename std::conditional<(RANK==0),
         typename std::conditional<(Access==Intent::ReadOnly),
             const value_type&,
             value_type&
         >::type,
         LocalView<value_type,RANK>
         >::type;
    };

    ArraySlicer(ArrayView<Value, Rank, AccessMode>& av) :
      av_(av) {
    }

    template <typename ...Args>
    typename Slice< SliceRank<Rank,Args...>::value >::type apply(const Args ...args) const {
        using ReturnType = typename Slice< SliceRank<Rank,Args...>::value >::type;
        return Slicer< ReturnType, (SliceRank<Rank,Args...>::value==0) >::apply( av_, args... );
    }

private:

    template <typename ReturnType, bool ToScalar = false>
    struct Slicer {
        template< typename ...Args >
        static ReturnType apply(ArrayView<Value, Rank, AccessMode>& av, const Args...args) {
            return ReturnType(
                  av.data()+offset(av,args...),
                  shape(av,args...).data(),
                  strides(av,args...).data()
            );
        }
    };

    template <typename ReturnType>
    struct Slicer<ReturnType,true> {
      template< typename ...Args >
      static ReturnType apply(ArrayView<Value, Rank, AccessMode>& av, const Args...args) {
          return *(av.data()+offset(av,args...));
      }
    };

    template< typename Int >
    static int start(Int idx) {
      return idx;
    }

    static int start(Range range) {
      return range.start;
    }

    template < int Dim, typename Int, typename... Ints >
    static constexpr int offset_part(ArrayView<Value, Rank, AccessMode>& av, const Int idx, const Ints... next_idx) {
        return start(idx)*av.stride(Dim) + offset_part<Dim+1>( av, next_idx... );
    }

    template < int Dim, typename Int >
    static constexpr int offset_part(ArrayView<Value, Rank, AccessMode>& av, const Int last_idx) {
        return start(last_idx)*av.stride(Dim);
    }

    template < typename... Args >
    static constexpr int offset(ArrayView<Value, Rank, AccessMode>& av, const Args... args) {
        return offset_part<0>(av,args...);
    }


    template< typename Shape, typename Int >
    static void update_shape( Shape& shape, int& i, const Int& ){
      // do nothing
    }
    template< typename Shape >
    static void update_shape( Shape& shape, int& i, const Range range ){
       shape[i++] = range.end-range.start;
    }

    template< int Dim, typename Shape, typename Int, typename... Ints >
    static void shape_part( ArrayView<Value, Rank, AccessMode>& av, Shape& shape, int& i, const Int idx, const Ints... next_idx ) {
      update_shape(shape,i,idx);
      shape_part<Dim+1>(av,shape,i,next_idx...);
    }

    template< int Dim, typename Shape, typename Int >
    static void shape_part( ArrayView<Value, Rank, AccessMode>& , Shape& shape, int& i, const Int idx) {
      update_shape(shape,i,idx);
    }

    template< typename... Args >
    static std::array<size_t,SliceRank<Rank,Args...>::value> shape(ArrayView<Value, Rank, AccessMode>& av, const Args... args ) {
      std::array<size_t,SliceRank<Rank,Args...>::value> result;
      int i(0);
      shape_part<0>(av,result,i,args...);
      return result;
    }


    template<int Dim, typename Strides, typename Int >
    static void update_strides( ArrayView<Value, Rank, AccessMode>&, Strides&, int& /*i*/, const Int& /*idx*/) {
      // do nothing
     }
    template<int Dim,typename Strides >
    static void update_strides( ArrayView<Value, Rank, AccessMode>& av, Strides& strides, int& i, const Range& /*range*/ ) {
       strides[i++] = av.stride(Dim);
    }

    template< int Dim, typename Strides, typename Int, typename... Ints >
    static void strides_part( ArrayView<Value, Rank, AccessMode>& av, Strides& strides, int& i, const Int idx, const Ints... next_idx ) {
      update_strides<Dim>(av,strides,i,idx);
      strides_part<Dim+1>(av,strides,i,next_idx...);
    }

    template< int Dim, typename Strides, typename Int >
    static void strides_part( ArrayView<Value, Rank, AccessMode>& av, Strides& strides, int& i, const Int idx) {
      update_strides<Dim>(av,strides,i,idx);
    }

    template< typename... Args >
    static std::array<size_t,SliceRank<Rank,Args...>::value> strides(ArrayView<Value, Rank, AccessMode>& av, const Args... args ) {
      std::array<size_t,SliceRank<Rank,Args...>::value> result;
      int i(0);
      strides_part<0>(av,result,i,args...);
      return result;
    }

private:
    ArrayView<Value, Rank, AccessMode>& av_;
};

//template <typename Value, int Rank, Intent AccessMode, typename ReturnType >
//struct ArraySlicer<Value, Rank, AccessMode, ReturnType, true> {
//    ArraySlicer(ArrayView<Value, Rank, AccessMode> const& av) : av_(av) {}
//    ArrayView<Value, Rank, AccessMode> const& av_;
//    ReturnType apply(const size_t i) const {
//        return *(av_.data_ + av_.strides_[0] * i);
//    }
//};
//------------------------------------------------------------------------------

} // namespace helpers
} // namespace array
} // namespace atlas
