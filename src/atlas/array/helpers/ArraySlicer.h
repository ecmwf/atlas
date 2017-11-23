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
#include "atlas/array.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace helpers {

class RangeBase {};

class From : public RangeBase {
public:
  From(int start) : start_(start) {}

  int start() const { return start_; }

  template < int Dim, typename View >
  int end(const View& view) const { return view.shape(Dim); }

private:
  int start_;
};

class To : public RangeBase {
public:
  To(int end) : end_(end) {}

  int start() const { return 0; }
  int end() const { return end_; }

private:
  int end_;
};

class All : public RangeBase {
public:
  All() {}

  int start() const { return 0; }

  template < int Dim, typename View >
  int end(const View& view) const { return view.shape(Dim); }
};

class Range : public RangeBase{
public:
  Range(int start, int end ) : start_(start), end_(end) {}
  int start() const { return start_; }
  int end() const { return end_; }

  static From from(int start) { return From(start); }
  static To   to(int end) { return To(end); }
  static All  all() { return All(); }

private:
  int start_;
  int end_;
};



class Index {
public:
  int index;
};

template< int Rank, typename ...Args >
struct SliceRank_impl;




template< int Dim >
struct deduce_slice_rank;


template<>
struct deduce_slice_rank<1> {

  template<typename Last>
  static constexpr int apply() {
    return std::is_base_of<RangeBase,Last>::value;
  }

};

#define EXPLICIT_TEMPLATE_SPECIALISATION(RANK) \
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
struct SliceRank_impl<1,Args...>
{
  enum { value = deduce_slice_rank<1>::apply<Args...>() };
};

template<typename Value, int RANK, Intent Access>
struct get_slice_type {
  using type =
    typename std::conditional<(RANK==0),
       typename std::conditional<(Access==Intent::ReadOnly),
           const Value&,
           Value&
       >::type,
       LocalView<Value,RANK>
    >::type;
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
      using type = typename get_slice_type< typename View::value_type, SliceRank<Args...>::value, View::ACCESS >::type;
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
          return *(view.data()+offset(view,args...));
      }
    };

    template< typename Int >
    static int start( Int idx ) {
      return idx;
    }

    static int start( Range range ) {
      return range.start();
    }

    static int start( All range ) {
      return range.start();
    }

    static int start( To range ) {
      return range.start();
    }

    static int start( From range ) {
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
    static void update_shape( View& view, Shape& shape, int& i, const All range ){
       shape[i++] = range.end<Dim>(view)-range.start();
    }
    template< int Dim, typename Shape >
    static void update_shape( View& view, Shape& shape, int& i, const From range ){
       shape[i++] = range.end<Dim>(view)-range.start();
    }
    template< int Dim, typename Shape >
    static void update_shape( View&, Shape& shape, int& i, const To range ){
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
    static void update_strides( View& view, Strides& strides, int& i, const From& /*range*/ ) {
       strides[i++] = view.stride(Dim);
    }
    template<int Dim,typename Strides >
    static void update_strides( View& view, Strides& strides, int& i, const To& /*range*/ ) {
       strides[i++] = view.stride(Dim);
    }
    template<int Dim,typename Strides >
    static void update_strides( View& view, Strides& strides, int& i, const All& /*range*/ ) {
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
