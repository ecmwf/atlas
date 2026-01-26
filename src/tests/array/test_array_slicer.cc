/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/helpers/ArraySlicer.h"
#include "atlas/library/config.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;
using namespace atlas::array::helpers;

template <typename Value, int Rank>
struct Slice {
    using type = typename std::conditional<(Rank == 0), Value, std::array<Value, Rank>>::type;
};

struct Slicer {
    template <int Rank, typename... Args>
    static typename Slice<double, SliceRank_impl<Rank, Args...>::value>::type apply(Args... args) {
        typename Slice<double, SliceRank_impl<Rank, Args...>::value>::type a;
    }
};

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_SliceRank") {
    static_assert(SliceRank_impl<1, int>::value == 0, "failed");
    static_assert(SliceRank_impl<1, Range>::value == 1, "failed");

    static_assert(SliceRank_impl<2, int, int>::value == 0, "failed");
    static_assert(SliceRank_impl<2, Range, int>::value == 1, "failed");
    static_assert(SliceRank_impl<2, Range, Range>::value == 2, "failed");

    static_assert(SliceRank_impl<3, int, int, int>::value == 0, "failed");
    static_assert(SliceRank_impl<3, Range, int, int>::value == 1, "failed");
    static_assert(SliceRank_impl<3, Range, Range, int>::value == 2, "failed");
    static_assert(SliceRank_impl<3, int, int, Range>::value == 1, "failed");
    static_assert(SliceRank_impl<3, Range, int, Range>::value == 2, "failed");
    static_assert(SliceRank_impl<3, Range, Range, Range>::value == 3, "failed");
    static_assert(SliceRank_impl<3, Range, int, Range>::value == 2, "failed");
    static_assert(SliceRank_impl<3, int, int, Range>::value == 1, "failed");
    static_assert(SliceRank_impl<3, int, Range, Range>::value == 2, "failed");
}

#if 1
CASE("test_array_slicer_1d") {
    ArrayT<double> arr(10);
    auto view = make_view<double, 1>(arr);
    view.assign({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

    using View = decltype(view);
    ArraySlicer<View> slicer(view);

    {
        auto slice = slicer.apply(Range{0, 2});
        static_assert(std::is_same<decltype(slice), LocalView<double, 1>>::value, "failed");
        EXPECT(slice.rank() == 1);
        EXPECT(slice.shape(0) == 2);
        EXPECT(slice(0) == 0);
        EXPECT(slice(1) == 1);
    }

    {
        auto slice = slicer.apply(Range{5, 10});
        static_assert(std::is_same<decltype(slice), LocalView<double, 1>>::value, "failed");
        EXPECT(slice.rank() == 1);
        EXPECT(slice.shape(0) == 5);
        EXPECT(slice(0) == 5);
        EXPECT(slice(1) == 6);
        EXPECT(slice(2) == 7);
        EXPECT(slice(3) == 8);
        EXPECT(slice(4) == 9);
    }

    {
        auto slice = slicer.apply(5);
        EXPECT(slice == 5);
        static_assert(std::is_same<decltype(slice), Reference<double>>::value, "failed");
    }

    {
        auto slice = slicer.apply(Range::all(), Range::dummy());
        static_assert(std::is_same<decltype(slice), LocalView<double, 2>>::value, "failed");
        EXPECT(slice.rank() == 2);
        EXPECT(slice.shape(0) == 10);
        EXPECT(slice.shape(1) == 1);
        EXPECT(slice(0, 0) == 0);
        EXPECT(slice(1, 0) == 1);
    }
}
#endif
CASE("test_array_slicer_2d") {
    ArrayT<double> arr(3, 5);
    auto view = make_view<double, 2>(arr);
    view.assign({
        11,
        12,
        13,
        14,
        15,
        21,
        22,
        23,
        24,
        25,
        31,
        32,
        33,
        34,
        35,
    });
    using View = decltype(view);
    ArraySlicer<View> slicer(view);

    //  ArraySlicer<double,2,Intent::ReadWrite> slicer(view);

    {
        auto slice = slicer.apply(Range{0, 2}, 2);
        EXPECT(slice.rank() == 1);
        EXPECT(slice.shape(0) == 2);
        EXPECT(slice(0) == 13);
        EXPECT(slice(1) == 23);
    }

    {
        const double& slice = slicer.apply(1, 3);
        EXPECT(slice == 24);
        // static_assert( std::is_same< decltype(slice) , double& >::value, "failed"
        // );
    }

    {
        auto slice = slicer.apply(Range{0, 2}, Range::dummy(), Range{1, 3});
        static_assert(std::is_same<decltype(slice), LocalView<double, 3>>::value, "failed");
        EXPECT(slice.rank() == 3);
        EXPECT(slice.shape(0) == 2);
        EXPECT(slice.shape(1) == 1);
        EXPECT(slice.shape(2) == 2);
        EXPECT(slice(0, 0, 0) == 12);
        EXPECT(slice(0, 0, 1) == 13);
        EXPECT(slice(1, 0, 0) == 22);
        EXPECT(slice(1, 0, 1) == 23);
    }
}

CASE("test_array_slicer_3d") {
    ArrayT<double> arr(2, 3, 5);
    auto view = make_view<double, 3>(arr);
    view.assign({111, 112, 113, 114, 115, 121, 122, 123, 124, 125, 131, 132, 133, 134, 135,

                 211, 212, 213, 214, 215, 221, 222, 223, 224, 225, 231, 232, 233, 234, 235});

    using View = decltype(view);
    ArraySlicer<View> slicer(view);

    {
        auto slice = slicer.apply(Range{0, 2}, 2, Range{2, 5});
        EXPECT(slice.rank() == 2);
        EXPECT(slice.shape(0) == 2);
        EXPECT(slice.shape(1) == 3);
        EXPECT(slice(0, 0) == 133);
        EXPECT(slice(0, 1) == 134);
        EXPECT(slice(0, 2) == 135);
        EXPECT(slice(1, 0) == 233);
        EXPECT(slice(1, 1) == 234);
        EXPECT(slice(1, 2) == 235);
    }

    {
        auto slice = slicer.apply(0, 0, Range::from(3));
        EXPECT(slice.rank() == 1);
        EXPECT(slice.shape(0) == 2);
        EXPECT(slice(0) == 114);
        EXPECT(slice(1) == 115);
    }

    {
        auto slice = slicer.apply(1, Range::all(), 3);
        EXPECT(slice.rank() == 1);
        EXPECT(slice.shape(0) == 3);
        EXPECT(slice(0) == 214);
        EXPECT(slice(1) == 224);
        EXPECT(slice(2) == 234);
    }

    {
        const double& slice = slicer.apply(1, 2, 3);
        EXPECT(slice == 234);
    }
}

CASE("test_array_slicer_of_slice") {
    ArrayT<double> arr(2, 3, 5);
    auto view = make_view<double, 3>(arr);
    view.assign({111, 112, 113, 114, 115, 121, 122, 123, 124, 125, 131, 132, 133, 134, 135,

                 211, 212, 213, 214, 215, 221, 222, 223, 224, 225, 231, 232, 233, 234, 235});

    using View = decltype(view);
    ArraySlicer<View> slicer(view);

    auto slice = slicer.apply(Range{0, 2}, 2, Range{2, 5});

    using Slice = decltype(slice);
    ArraySlicer<Slice> subslicer(slice);

    auto subslice1 = subslicer.apply(0, 0);
    auto subslice2 = subslicer.apply(Range::all(), Range::all());
    auto subslice3 = subslicer.apply(Range::all(), 0);
    auto subslice4 = subslicer.apply(0, Range::to(2));

    EXPECT(is_approximately_equal(double(subslice1), 133.));
    EXPECT(subslice2.size() == slice.size());
    EXPECT(subslice3.size() == slice.shape(0));
    EXPECT(subslice4.size() == 2);
}

CASE("test_arrayview_slice_type") {
    ArrayT<double> arr(2, 3, 5);

    // Assign
    {
        auto view = make_view<double, 3>(arr);
        view.assign({111, 112, 113, 114, 115, 121, 122, 123, 124, 125, 131, 132, 133, 134, 135,

                     211, 212, 213, 214, 215, 221, 222, 223, 224, 225, 231, 232, 233, 234, 235});
    }

    {
        auto read_write_view = make_view<double, 3>(arr);

        auto slice1 = read_write_view.slice(Range{0, 2}, 2, Range{2, 5});

        EXPECT(slice1(0, 0) == 133);
        EXPECT(slice1(0, 1) == 134);
        EXPECT(slice1(0, 2) == 135);
        EXPECT(slice1(1, 0) == 233);
        EXPECT(slice1(1, 1) == 234);
        EXPECT(slice1(1, 2) == 235);

        auto slice2 = read_write_view.slice(Range::all(), Range::to(2), Range::from(3));

        EXPECT(slice2(0, 0, 0) == 114);
        EXPECT(slice2(0, 0, 1) == 115);
        EXPECT(slice2(0, 1, 0) == 124);
        EXPECT(slice2(0, 1, 1) == 125);
        EXPECT(slice2(1, 0, 0) == 214);
        EXPECT(slice2(1, 0, 1) == 215);
        EXPECT(slice2(1, 1, 0) == 224);
        EXPECT(slice2(1, 1, 1) == 225);

        auto slice3 = read_write_view.slice(0, 1, 2);

        EXPECT(slice3 == 123);

        static_assert(std::is_same<decltype(slice1), LocalView<double, 2>>::value, "failed");
        static_assert(std::is_same<decltype(slice2), LocalView<double, 3>>::value, "failed");
        static_assert(std::is_same<decltype(slice3), Reference<double>>::value, "failed");
    }

    {
        auto read_only_view = make_view<const double, 3>(arr);

        auto slice1 = read_only_view.slice(Range{0, 2}, 2, Range{2, 5});

        EXPECT(slice1(0, 0) == 133);
        EXPECT(slice1(0, 1) == 134);
        EXPECT(slice1(0, 2) == 135);
        EXPECT(slice1(1, 0) == 233);
        EXPECT(slice1(1, 1) == 234);
        EXPECT(slice1(1, 2) == 235);

        auto slice2 = read_only_view.slice(Range::all(), Range::to(2), Range::from(3));

        EXPECT(slice2(0, 0, 0) == 114);
        EXPECT(slice2(0, 0, 1) == 115);
        EXPECT(slice2(0, 1, 0) == 124);
        EXPECT(slice2(0, 1, 1) == 125);
        EXPECT(slice2(1, 0, 0) == 214);
        EXPECT(slice2(1, 0, 1) == 215);
        EXPECT(slice2(1, 1, 0) == 224);
        EXPECT(slice2(1, 1, 1) == 225);

        auto slice3 = read_only_view.slice(0, 1, 2);

        EXPECT(slice3 == 123);

        const double& slice4 = read_only_view.slice(0, 1, 2);

        EXPECT(slice4 == 123);

        static_assert(std::is_same<decltype(slice1), LocalView<const double, 2>>::value, "failed");
        static_assert(std::is_same<decltype(slice2), LocalView<const double, 3>>::value, "failed");
        static_assert(std::is_same<decltype(slice3), Reference<const double>>::value, "failed");
        static_assert(std::is_same<decltype(slice4), double const&>::value, "failed");
    }

    {
        const auto const_read_write_view = make_view<double, 3>(arr);

        auto slice1 = const_read_write_view.slice(Range{0, 2}, 2, Range{2, 5});

        EXPECT(slice1(0, 0) == 133);
        EXPECT(slice1(0, 1) == 134);
        EXPECT(slice1(0, 2) == 135);
        EXPECT(slice1(1, 0) == 233);
        EXPECT(slice1(1, 1) == 234);
        EXPECT(slice1(1, 2) == 235);

        auto slice2 = const_read_write_view.slice(Range::all(), Range::to(2), Range::from(3));

        EXPECT(slice2(0, 0, 0) == 114);
        EXPECT(slice2(0, 0, 1) == 115);
        EXPECT(slice2(0, 1, 0) == 124);
        EXPECT(slice2(0, 1, 1) == 125);
        EXPECT(slice2(1, 0, 0) == 214);
        EXPECT(slice2(1, 0, 1) == 215);
        EXPECT(slice2(1, 1, 0) == 224);
        EXPECT(slice2(1, 1, 1) == 225);

        auto slice3 = const_read_write_view.slice(0, 1, 2);

        EXPECT(slice3 == 123);

        static_assert(std::is_same<decltype(slice1), LocalView<const double, 2>>::value, "failed");
        static_assert(std::is_same<decltype(slice2), LocalView<const double, 3>>::value, "failed");
        static_assert(std::is_same<decltype(slice3), Reference<double const>>::value, "failed");
    }


    {
        const auto const_read_write_view = make_view<double, 3>(arr);

        auto slice1 = const_read_write_view.slice(Range{0, 2}, 2, Range{2, 5}).as_mdspan();
#define INDEX(...) std::array{__VA_ARGS__}
        EXPECT(slice1[INDEX(0, 0)] == 133);
        EXPECT(slice1[INDEX(0, 1)] == 134);
        EXPECT(slice1[INDEX(0, 2)] == 135);
        EXPECT(slice1[INDEX(1, 0)] == 233);
        EXPECT(slice1[INDEX(1, 1)] == 234);
        EXPECT(slice1[INDEX(1, 2)] == 235);

        auto slice2 = const_read_write_view.slice(Range::all(), Range::to(2), Range::from(3)).as_mdspan();

        EXPECT(slice2[INDEX(0, 0, 0)] == 114);
        EXPECT(slice2[INDEX(0, 0, 1)] == 115);
        EXPECT(slice2[INDEX(0, 1, 0)] == 124);
        EXPECT(slice2[INDEX(0, 1, 1)] == 125);
        EXPECT(slice2[INDEX(1, 0, 0)] == 214);
        EXPECT(slice2[INDEX(1, 0, 1)] == 215);
        EXPECT(slice2[INDEX(1, 1, 0)] == 224);
        EXPECT(slice2[INDEX(1, 1, 1)] == 225);

        auto slice2_view = LocalView<const double,3>(slice2);

        EXPECT(slice2_view(0, 0, 0) == 114);
        EXPECT(slice2_view(0, 0, 1) == 115);
        EXPECT(slice2_view(0, 1, 0) == 124);
        EXPECT(slice2_view(0, 1, 1) == 125);
        EXPECT(slice2_view(1, 0, 0) == 214);
        EXPECT(slice2_view(1, 0, 1) == 215);
        EXPECT(slice2_view(1, 1, 0) == 224);
        EXPECT(slice2_view(1, 1, 1) == 225);

        static_assert(std::is_same<decltype(slice1), mdspan<const double, dims<2>, layout_stride>>::value, "failed");
        static_assert(std::is_same<decltype(slice2), mdspan<const double, dims<3>, layout_stride>>::value, "failed");
    }


}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
