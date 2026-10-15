#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "atlas/array/Array.h"
#include "atlas/array/LocalView.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/native/NativeIndexView.h"
#include "atlas/array/mdspan.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;

using namespace atlas::array::introspection;

namespace atlas {
namespace test {

using DynamicFirstFixedLastExtents = extents<std::size_t, dynamic_extent, 3>;

template <typename Span, typename View>
void expect_same_data(const Span& span, const View& view) {
    static_assert(std::is_same_v<typename Span::accessor_type, default_accessor<typename Span::element_type>>);
    using DefaultSpan = mdspan<typename Span::element_type, typename Span::extents_type, typename Span::layout_type>;
    static_assert(std::is_same_v<typename Span::accessor_type, typename DefaultSpan::accessor_type>);
    EXPECT(span.data_handle() == view.data());
}

template <typename Span, typename View>
void expect_shape_2d(const Span& span, const View& view) {
    EXPECT(span.extent(0) == view.shape(0));
    EXPECT(span.extent(1) == view.shape(1));
}

template <typename Span, typename View>
void expect_strides_2d(const Span& span, const View& view) {
    EXPECT(span.stride(0) == view.stride(0));
    EXPECT(span.stride(1) == view.stride(1));
}

template <typename Span>
void expect_dynamic_first_fixed_last_extents() {
    static_assert(Span::extents_type::static_extent(0) == dynamic_extent);
    static_assert(Span::extents_type::static_extent(1) == 3);
}

template <typename Span, typename View>
void expect_layout_stride_2d(const Span& span, const View& view) {
    static_assert(std::is_same_v<typename Span::layout_type, layout_stride>);
    expect_same_data(span, view);
    expect_shape_2d(span, view);
    expect_strides_2d(span, view);
}

template <typename Span, typename View>
void expect_layout_right_2d(const Span& span, const View& view) {
    static_assert(std::is_same_v<typename Span::layout_type, layout_right>);
    expect_same_data(span, view);
    expect_shape_2d(span, view);
}

LocalView<double, 2> make_local_view(double (&data)[6]) {
    idx_t shape[2]{2, 3};
    idx_t strides[2]{3, 1};
    return LocalView<double, 2>{data, shape, strides};
}

LocalView<double, 2> make_local_view_4x3(double (&data)[12]) {
    idx_t shape[2]{4, 3};
    idx_t strides[2]{3, 1};
    return LocalView<double, 2>{data, shape, strides};
}

struct LayoutRightViewWithoutStrides {
    using layout_type = layout_right;
    static constexpr std::size_t rank() { return 2; }
};

struct OffsetAccessor : default_accessor<double> {
    std::size_t displacement{0};

    double& access(double* data, std::size_t offset) const noexcept {
        return data[offset + displacement];
    }
};

struct RequiredOffsetAccessor : OffsetAccessor {
    RequiredOffsetAccessor() = delete;
    explicit RequiredOffsetAccessor(std::size_t offset) { displacement = offset; }
};

template <typename Accessor>
struct ViewWithAccessor {
    using element_type = double;
    using data_handle_type = double*;
    using accessor_type = Accessor;

    double* data;
    Accessor accessor_instance;

    static constexpr std::size_t rank() { return 1; }
    std::size_t extent(std::size_t) const { return 2; }
    std::size_t stride(std::size_t) const { return 1; }
    double* data_handle() const { return data; }
    const Accessor& accessor() const { return accessor_instance; }
};

CASE("test_make_mdspan_custom_view_preserves_accessor") {
    double data[4]{10., 20., 30., 40.};
    auto check = [&](const auto& view) {
        using View = std::decay_t<decltype(view)>;
        static_assert(!introspection::is_mdspan<View>());
        auto check_span = [&](auto span) {
            static_assert(std::is_same_v<typename decltype(span)::accessor_type, typename View::accessor_type>);
            EXPECT(span.data_handle() == data);
            EXPECT(span.extent(0) == 2);
            EXPECT(span.stride(0) == 1);
            EXPECT(span.accessor().displacement == 1);
            EXPECT(span(0) == 20.);
            EXPECT(span(1) == 30.);
        };
        check_span(make_mdspan(view));
        check_span(make_mdspan<layout_right>(view));
        check_span(make_mdspan(view, std::array<std::size_t, 1>{2}));
        check_span(make_mdspan<dims<1>>(view));
        check_span(make_mdspan<dims<1>>(view, std::array<std::size_t, 1>{2}));
        auto replaced = make_mdspan<default_accessor>(view);
        static_assert(std::is_same_v<typename decltype(replaced)::accessor_type, default_accessor<double>>);
        EXPECT(replaced(0) == 10.);
    };
    OffsetAccessor accessor;
    accessor.displacement = 1;
    check(ViewWithAccessor<OffsetAccessor>{data, accessor});
    check(ViewWithAccessor<RequiredOffsetAccessor>{data, RequiredOffsetAccessor{1}});
}

CASE("test_make_mdspan_invalid_argument_details") {
    double data[8]{};
    idx_t shape[2]{2, 3};
    auto check = [](auto operation, std::initializer_list<const char*> details) {
        bool caught = false;
        try {
            operation();
        }
        catch (const std::invalid_argument& error) {
            caught = true;
            const std::string message = error.what();
            for (const auto* detail : details) {
                EXPECT(message.find(detail) != std::string::npos);
            }
        }
        EXPECT(caught);
    };
    idx_t padded_strides[2]{4, 1};
    LocalView<double, 2> padded{data, shape, padded_strides};
    check([&] { (void)make_mdspan<layout_right>(padded); },
          {"source strides are incompatible with layout_right",
           "extents={2, 3}", "source strides={4, 1}", "required strides={3, 1}"});
    check([&] { (void)make_mdspan<layout_left>(padded); },
          {"source strides are incompatible with layout_left",
           "extents={2, 3}", "source strides={4, 1}", "required strides={1, 2}"});
    idx_t negative_strides[2]{-3, 1};
    LocalView<double, 2> negative{data, shape, negative_strides};
    check([&] { (void)make_mdspan(negative); }, {"negative source stride", "dimension 0", "-3"});
    idx_t zero_strides[2]{0, 1};
    LocalView<double, 2> zero{data, shape, zero_strides};
    check([&] { (void)make_mdspan(zero); },
          {"layout_stride", "positive, non-overlapping strides", "extents={2, 3}", "strides={0, 1}"});
    idx_t overlapping_strides[2]{1, 1};
    LocalView<double, 2> overlapping{data, shape, overlapping_strides};
    check([&] { (void)make_mdspan(overlapping); },
          {"layout_stride", "positive, non-overlapping strides", "extents={2, 3}", "strides={1, 1}"});
}

CASE("test_make_mdspan_layout_compatibility") {
    double data[8]{10., 20., 30., 40., 50., 60., 70., 80.};
    using Extents = extents<std::size_t, 2, 3>;
    using LeftSpan = mdspan<double, Extents, layout_left, OffsetAccessor>;
    OffsetAccessor accessor;
    accessor.displacement = 1;
    LeftSpan source{data, typename LeftSpan::mapping_type{}, accessor};
    auto preserved = make_mdspan(source);
    auto shaped = make_mdspan(source, std::array<std::size_t, 2>{2, 3});
    static_assert(std::is_same_v<decltype(preserved), LeftSpan>);
    EXPECT(preserved.data_handle() == data);
    EXPECT(preserved.accessor().displacement == 1);
    EXPECT(preserved(1, 2) == source(1, 2));
    EXPECT(shaped(1, 2) == source(1, 2));
    auto strided = make_mdspan<layout_stride>(source);
    EXPECT(conforms_layout<layout_left>(strided));
    EXPECT(!conforms_layout<layout_right>(strided));
    auto restored = make_mdspan<layout_left>(strided);
    EXPECT(restored(1, 2) == source(1, 2));

    mdspan<double, Extents, layout_right> right{data};
    auto incompatible_shape = [&] { return make_mdspan<layout_right>(right, std::array<std::size_t, 2>{3, 2}); };
    EXPECT_THROWS_AS(incompatible_shape(), std::invalid_argument);
    EXPECT(!conforms_layout<layout_left>(right));
    EXPECT(!conforms_layout<layout_right>(source));

    idx_t shape[2]{2, 3};
    idx_t padded[2]{4, 1};
    LocalView<double, 2> padded_view{data, shape, padded};
    auto padded_right = [&] { return make_mdspan<layout_right>(padded_view); };
    auto padded_left = [&] { return make_mdspan<layout_left>(padded_view); };
    EXPECT_THROWS_AS(padded_right(), std::invalid_argument);
    EXPECT_THROWS_AS(padded_left(), std::invalid_argument);
    idx_t overlapping_strides[2]{1, 1};
    LocalView<double, 2> overlapping{data, shape, overlapping_strides};
    EXPECT(!conforms_layout<layout_stride>(overlapping));
    auto overlapping_conversion = [&] { return make_mdspan(overlapping); };
    EXPECT_THROWS_AS(overlapping_conversion(), std::invalid_argument);
    auto strided_right = make_mdspan<layout_stride>(right);
    auto wrong_order = [&] { return make_mdspan<layout_left>(strided_right); };
    EXPECT_THROWS_AS(wrong_order(), std::invalid_argument);

    idx_t singleton_shape[2]{1, 3};
    idx_t singleton_strides[2]{99, 1};
    LocalView<double, 2> singleton{data, singleton_shape, singleton_strides};
    EXPECT(conforms_layout<layout_left>(singleton));
    EXPECT(conforms_layout<layout_right>(singleton));
    EXPECT(make_mdspan<layout_left>(singleton)(0, 2) == data[2]);
    idx_t empty_shape[2]{2, 0};
    LocalView<double, 2> empty{data, empty_shape, padded};
    EXPECT(conforms_layout<layout_left>(empty));
    EXPECT(make_mdspan<layout_left>(empty).size() == 0);

    std::vector<double> values{1., 2., 3., 4., 5., 6.};
    auto check_container = [](auto& container) {
        auto left = make_mdspan<layout_left>(container, std::array<std::size_t, 2>{2, 3});
        static_assert(std::is_same_v<typename decltype(left)::layout_type, layout_left>);
        EXPECT(left.data_handle() == container.data());
        EXPECT(left.stride(0) == 1);
        EXPECT(left.stride(1) == 2);
        EXPECT(left(1, 0) == 2.);
        EXPECT(left(0, 2) == 5.);
    };
    check_container(values);
    check_container(std::as_const(values));
    std::array<double, 6> array_values{1., 2., 3., 4., 5., 6.};
    check_container(array_values);
    auto pointer = make_mdspan<layout_left>(values.data(), std::array<std::size_t, 2>{2, 3});
    EXPECT(pointer(1, 0) == 2.);
    std::vector<double> empty_values;
    auto empty_left = make_mdspan<layout_left>(empty_values, std::array<std::size_t, 2>{0, 3});
    EXPECT(empty_left.size() == 0);
    EXPECT(empty_left.extent(0) == 0);
    EXPECT(empty_left.extent(1) == 3);
}

CASE("test_make_mdspan_review_rank_one_strided_layout") {
    // A rank-1 view with stride two must not qualify for contiguous row-major conversion.
    double data[5]{1., 2., 3., 4., 5.};
    idx_t shape[1]{3};
    idx_t strides[1]{2};
    auto view = LocalView<double, 1>{data, shape, strides};
    auto span = make_mdspan<layout_stride>(view);
    EXPECT(span(1) == 3.);
    EXPECT(not conforms_layout_right(view));
    EXPECT(not conforms_layout_right(span));
    EXPECT(not conforms_layout<layout_right>(view));
    auto incompatible_right = [&] { return make_mdspan<layout_right>(view); };
    EXPECT_THROWS_AS(incompatible_right(), std::invalid_argument);

    idx_t contiguous_strides[1]{1};
    auto contiguous_view = LocalView<double, 1>{data, shape, contiguous_strides};
    EXPECT(conforms_layout_right(contiguous_view));
    auto contiguous_span = make_mdspan<layout_right>(contiguous_view);
    EXPECT(contiguous_span(1) == 2.);
    std::vector<double> values{1., 2., 3.};
    EXPECT(conforms_layout_right(values));
    EXPECT(conforms_layout_right(std::as_const(values)));
    std::array<double, 3> array_values{1., 2., 3.};
    EXPECT(conforms_layout_right(array_values));
}

CASE("test_make_mdspan_review_preserve_accessor_state") {
    // Conversion must retain the source accessor instance, not only its type.
    double data[4]{10., 20., 30., 40.};
    using Span = mdspan<double, extents<std::size_t, 2>, layout_right, OffsetAccessor>;
    OffsetAccessor accessor;
    accessor.displacement = 1;
    Span source{data, typename Span::mapping_type{}, accessor};
    EXPECT(source(0) == 20.);
    auto converted = make_mdspan(source);
    auto strided = make_mdspan<layout_stride>(source);
    auto shaped = make_mdspan(source, typename Span::extents_type{});
    auto explicit_extents = make_mdspan<dims<1>>(source);
    auto explicit_shaped = make_mdspan<dims<1>, layout_right>(source, std::array<std::size_t, 1>{2});
    EXPECT(converted.accessor().displacement == 1);
    EXPECT(converted(0) == source(0));
    EXPECT(strided(0) == source(0));
    EXPECT(shaped(0) == source(0));
    EXPECT(explicit_extents.accessor().displacement == 1);
    EXPECT(explicit_extents(0) == source(0));
    EXPECT(explicit_shaped.accessor().displacement == 1);
    EXPECT(explicit_shaped(0) == source(0));
    auto replaced = make_mdspan<default_accessor>(source);
    EXPECT(replaced(0) == data[0]);
    auto explicit_replaced = make_mdspan<dims<1>, layout_right, default_accessor>(source);
    auto explicit_shaped_replaced = make_mdspan<dims<1>, layout_right, default_accessor>(source, std::array<std::size_t, 1>{2});
    EXPECT(explicit_replaced(0) == data[0]);
    EXPECT(explicit_shaped_replaced(0) == data[0]);

    using RequiredSpan = mdspan<double, extents<std::size_t, 2>, layout_right, RequiredOffsetAccessor>;
    RequiredSpan required{data, typename RequiredSpan::mapping_type{}, RequiredOffsetAccessor{1}};
    auto required_copy = make_mdspan(std::as_const(required));
    auto required_strided = make_mdspan<layout_stride>(required);
    auto required_shaped = make_mdspan(required, typename RequiredSpan::extents_type{});
    auto required_explicit = make_mdspan<dims<1>, layout_right>(required);
    auto required_explicit_shaped = make_mdspan<dims<1>>(required, std::array<std::size_t, 1>{2});
    EXPECT(required_copy(0) == 20.);
    EXPECT(required_strided(0) == 20.);
    EXPECT(required_shaped(0) == 20.);
    EXPECT(required_explicit(0) == 20.);
    EXPECT(required_explicit_shaped(0) == 20.);
}

template <typename Source, typename = void>
struct CanMakeMdspan : std::false_type {};

template <typename Source>
struct CanMakeMdspan<Source, std::void_t<decltype(make_mdspan(std::declval<Source>()))>> : std::true_type {};

template <typename Source, typename = void>
struct CanMakeShapedMdspan : std::false_type {};

template <typename Source>
struct CanMakeShapedMdspan<Source, std::void_t<decltype(
    make_mdspan(std::declval<Source>(), std::declval<std::array<std::size_t, 2>>()))>> : std::true_type {};

CASE("test_make_mdspan_temporary_view_lifetime") {
    // Temporary non-owning views copy metadata safely; temporary owning containers are rejected.
    static_assert(std::is_trivially_copyable_v<LocalView<double, 2>>);
    static_assert(!CanMakeMdspan<std::vector<double>&&>::value);
    static_assert(!CanMakeMdspan<const std::vector<double>&&>::value);
    static_assert(!CanMakeMdspan<std::array<double, 6>&&>::value);
    static_assert(!CanMakeShapedMdspan<std::vector<double>&&>::value);
    static_assert(CanMakeMdspan<std::vector<double>&>::value);
    static_assert(CanMakeMdspan<LocalView<double, 2>&&>::value);

    double data[6]{1., 2., 3., 4., 5., 6.};
    auto original = make_local_view(data);
    auto copied = original;
    EXPECT(copied.shape() != original.shape());
    EXPECT(copied.strides() != original.strides());
    EXPECT(copied(1, 2) == 6.);
    copied = original;
    EXPECT(copied.shape() != original.shape());
    EXPECT(copied.strides() != original.strides());
    auto span = make_mdspan(make_local_view(data));
    EXPECT(span.data_handle() == data);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.stride(0) == 3);
    span(1, 2) = 12.;
    EXPECT(data[5] == 12.);

    auto converted = make_mdspan<layout_right>(make_mdspan(make_local_view(data)));
    EXPECT(converted(1, 2) == 12.);
    auto shaped = make_mdspan(make_local_view(data), std::array<std::size_t, 2>{2, 3});
    EXPECT(shaped(1, 2) == 12.);
    auto default_span = make_mdspan<default_accessor>(make_local_view(data));
    auto aligned_span = make_mdspan<aligned_accessor>(make_local_view(data));
    auto layout_span = make_mdspan<layout_right, default_accessor>(make_local_view(data));
    auto aligned_layout = make_mdspan<fixed_last_dim<3>, aligned_accessor>(make_local_view(data));
    auto fixed_span = make_mdspan<fixed_last_dim<3>, layout_right, default_accessor>(make_local_view(data));
    auto fixed_aligned = make_mdspan<fixed_last_dim<3>, layout_right, aligned_accessor>(make_local_view(data));
    EXPECT(default_span(1, 2) == 12.);
    EXPECT(aligned_span(1, 2) == 12.);
    EXPECT(layout_span(1, 2) == 12.);
    EXPECT(aligned_layout(1, 2) == 12.);
    EXPECT(fixed_span(1, 2) == 12.);
    EXPECT(fixed_aligned(1, 2) == 12.);
}

CASE("test_make_mdspan_as_const_views") {
    // as_const() preserves storage, shape, strides, and index adjustment while preventing writes.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);
    view(1, 2) = 42.;
    const auto& const_view = view;
    auto readonly = make_mdspan(const_view.as_const());
    auto readonly_again = make_mdspan(view.as_const().as_const());
    static_assert(std::is_same_v<typename decltype(readonly)::element_type, const double>);
    static_assert(!std::is_assignable_v<decltype(readonly(1, 2)), double>);
    EXPECT(readonly.data_handle() == view.data());
    expect_shape_2d(readonly, view);
    expect_strides_2d(readonly, view);
    view(1, 2) = 43.;
    EXPECT(readonly(1, 2) == 43.);
    EXPECT(readonly_again(1, 2) == 43.);

    double data[8]{};
    idx_t shape[2]{2, 3};
    idx_t strides[2]{4, 1};
    auto local = LocalView<double, 2>{data, shape, strides};
    local(1, 2) = 17.;
    auto readonly_local = make_mdspan(local.as_const());
    auto readonly_local_again = make_mdspan(local.as_const().as_const());
    static_assert(std::is_same_v<typename decltype(readonly_local)::element_type, const double>);
    static_assert(!std::is_assignable_v<decltype(readonly_local(1, 2)), double>);
    EXPECT(readonly_local.data_handle() == data);
    expect_shape_2d(readonly_local, local);
    expect_strides_2d(readonly_local, local);
    EXPECT(readonly_local(1, 2) == 17.);
    EXPECT(readonly_local_again(1, 2) == 17.);

    int indices[8]{};
    idx_t index_shape[1]{4};
    idx_t index_strides[1]{2};
    auto index_view = IndexView<int, 1>{indices, index_shape, index_strides};
    index_view(3) = 7;
    auto readonly_index = make_mdspan(index_view.as_const());
    auto readonly_index_again = make_mdspan(index_view.as_const().as_const());
    static_assert(std::is_same_v<typename decltype(readonly_index)::element_type, const int>);
    static_assert(std::is_same_v<typename decltype(readonly_index)::accessor_type,
                                index_accessor<const int, ATLAS_HAVE_FORTRAN>>);
    static_assert(!std::is_assignable_v<decltype(readonly_index(3)), int>);
    EXPECT(readonly_index.data_handle() == indices);
    EXPECT(readonly_index.extent(0) == index_view.extent(0));
    EXPECT(readonly_index.stride(0) == index_view.stride(0));
    EXPECT(readonly_index(3) == 7);
    EXPECT(indices[6] == 7 + ATLAS_HAVE_FORTRAN);
    auto explicit_index = make_mdspan<dims<1>>(index_view);
    auto explicit_shaped_index = make_mdspan<dims<1>>(index_view, std::array<std::size_t, 1>{4});
    static_assert(std::is_same_v<typename decltype(explicit_index)::accessor_type,
                                index_accessor<int, ATLAS_HAVE_FORTRAN>>);
    static_assert(std::is_same_v<typename decltype(explicit_shaped_index)::accessor_type,
                                index_accessor<int, ATLAS_HAVE_FORTRAN>>);
    EXPECT(explicit_index(3) == 7);
    EXPECT(explicit_shaped_index(3) == 7);
    index_view(3) = 8;
    EXPECT(readonly_index_again(3) == 8);
}

CASE("test_make_mdspan_arrayview_layout_stride") {
    // ArrayView conversion preserves its data, shape, and strides with layout_stride.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 2) = 12.;

    auto span = make_mdspan<layout_stride>(view);

    expect_layout_stride_2d<decltype(span)>(span, view);
    EXPECT(span(1, 2) == 12.);
}

CASE("test_make_mdspan_arrayview_layout_right") {
    // A contiguous ArrayView can be exposed with layout_right without copying its elements.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 2) = 42.;

    auto span = make_mdspan<layout_right>(view);

    expect_layout_right_2d<decltype(span)>(span, view);
    EXPECT(span(1, 2) == 42.);
}

CASE("test_make_mdspan_arrayview_free_function") {
    // The free conversion function aliases ArrayView storage and reads the original values.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 1) = 7.;

    auto span = make_mdspan<layout_stride>(view);

    expect_same_data(span, view);
    EXPECT(span(1, 1) == 7.);
}

CASE("test_array_is_mdspan") {
    // Type-based and object-based introspection distinguish Atlas views from mdspans.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);
    auto span = make_mdspan<layout_stride>(view);

    static_assert(!array::is_mdspan<decltype(view)>());
    static_assert(array::is_mdspan<decltype(span)>());
    static_assert(!array::is_mdspan(view));
    static_assert(array::is_mdspan(span));
    EXPECT(span.data_handle() == view.data());
}

CASE("test_make_mdspan_const_arrayview") {
    // Object constness is shallow; only a const element type makes the converted span read-only.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);
    view(0, 2) = 5.;

    const auto& const_view = view;
    auto span = make_mdspan<layout_stride>(const_view);

    static_assert(std::is_same_v<typename decltype(span)::element_type, double>);
    expect_same_data(span, const_view);
    EXPECT(span(0, 2) == 5.);
    span(0, 2) = 6.;
    EXPECT(view(0, 2) == 6.);
    auto explicit_span = make_mdspan<layout_stride, default_accessor>(const_view);
    static_assert(std::is_same_v<typename decltype(explicit_span)::element_type, double>);
    explicit_span(0, 2) = 7.;
    EXPECT(view(0, 2) == 7.);

    auto readonly_view = make_host_view<const double, 2>(*array);
    const auto& const_readonly_view = readonly_view;
    auto readonly_span = make_mdspan(const_readonly_view);
    static_assert(std::is_same_v<typename decltype(readonly_span)::element_type, const double>);
    static_assert(!std::is_assignable_v<decltype(readonly_span(0, 2)), double>);
    EXPECT(readonly_span(0, 2) == 7.);
}

CASE("test_make_mdspan_runtime_extents") {
    // Explicit runtime extents preserve the ArrayView shape and stride mapping.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    std::array<std::size_t, 2> extents{2, 3};
    auto span = make_mdspan<layout_stride>(view, extents);

    expect_layout_stride_2d<decltype(span)>(span, view);
}

CASE("test_make_mdspan_arrayview_dynamic_fixed_extents") {
    // Both layouts support a dynamic outer extent and a compile-time fixed inner extent.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 2) = 15.;

    auto stride_span = make_mdspan<DynamicFirstFixedLastExtents, layout_stride>(view);
    auto right_span = make_mdspan<DynamicFirstFixedLastExtents, layout_right>(view);

    expect_dynamic_first_fixed_last_extents<decltype(stride_span)>();
    expect_layout_stride_2d<decltype(stride_span)>(stride_span, view);
    EXPECT(stride_span(1, 2) == 15.);
    expect_layout_right_2d<decltype(right_span)>(right_span, view);
    EXPECT(right_span(1, 2) == 15.);
}

CASE("test_make_mdspan_dynamic_fixed_extents_argument") {
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    // An extents object supplies the runtime outer size while retaining the fixed inner size.
    auto view = make_host_view<double, 2>(*array);

    view(1, 2) = 18.;

    DynamicFirstFixedLastExtents extents{2};
    auto stride_span = make_mdspan<layout_stride>(view, extents);
    auto right_span = make_mdspan<layout_right>(view, extents);

    expect_dynamic_first_fixed_last_extents<decltype(stride_span)>();
    expect_layout_stride_2d<decltype(stride_span)>(stride_span, view);
    EXPECT(stride_span(1, 2) == 18.);
    expect_layout_right_2d<decltype(right_span)>(right_span, view);
    EXPECT(right_span(1, 2) == 18.);
}

CASE("test_make_mdspan_static_extent_transform") {
    // Dimension transforms refine static extents while allowing layout and accessor selection.
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 2) = 63.;

    auto dim_span = make_mdspan<fixed_dim<1, 3>>(view);
    auto last_span = make_mdspan<fixed_last_dim<3>, layout_right, default_accessor>(view);

    static_assert(decltype(dim_span)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(dim_span)::extents_type::static_extent(1) == 3);
    static_assert(std::is_same_v<typename decltype(dim_span)::layout_type, layout_stride>);
    static_assert(std::is_same_v<typename decltype(dim_span)::accessor_type, default_accessor<double>>);
    static_assert(decltype(last_span)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(last_span)::extents_type::static_extent(1) == 3);
    static_assert(std::is_same_v<typename decltype(last_span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(last_span)::accessor_type, default_accessor<double>>);
    EXPECT(dim_span.extent(0) == 2);
    EXPECT(dim_span.extent(1) == 3);
    EXPECT(dim_span(1, 2) == 63.);
    EXPECT(last_span(1, 2) == 63.);
}

CASE("test_make_mdspan_static_extent_transform_argument") {
    // A fixed-last-dimension transform also works with an explicitly supplied runtime shape.
    std::unique_ptr<Array> array{Array::create<double>(4, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(3, 2) = 64.;

    std::array<std::size_t, 2> extents{4, 3};
    auto span = make_mdspan<fixed_last_dim<3>, layout_right, default_accessor>(view, extents);

    static_assert(decltype(span)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(span)::extents_type::static_extent(1) == 3);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, default_accessor<double>>);
    EXPECT(span.extent(0) == 4);
    EXPECT(span.extent(1) == 3);
    EXPECT(span(3, 2) == 64.);
}

CASE("test_array_static_extent_introspection") {
    // All static-extent query forms distinguish dynamic extents from compile-time fixed sizes.
    std::unique_ptr<Array> array{Array::create<double>(4, 3)};
    auto view = make_host_view<double, 2>(*array);

    auto span = make_mdspan<fixed_last_dim<3>, layout_right>(view);

    static_assert(array::static_extent<decltype(span)>(0) == dynamic_extent);
    static_assert(array::static_extent<decltype(span)>(1) == 3);
    static_assert(array::static_extent(span, 0) == dynamic_extent);
    static_assert(array::static_extent(span, 1) == 3);
    static_assert(array::static_extent<0>(span) == dynamic_extent);
    static_assert(array::static_extent<1>(span) == 3);
    EXPECT(span.data_handle() == view.data());
    EXPECT(span.extent(0) == 4);
    EXPECT(span.extent(1) == 3);
}

CASE("test_make_mdspan_static_last_extent_layout_right_restrict_aligned_policy") {
    // A blocked aligned view supports a fixed SIMD extent and a restrict-aligned accessor policy.
    constexpr idx_t nproma = 32;
    alignas(64) double data[2 * 4 * nproma]{};
    idx_t shape[3]{2, 4, nproma};
    idx_t strides[3]{4 * nproma, nproma, 1};
    auto view = LocalView<double, 3>{data, shape, strides};

    view(1, 3, 31) = 71.;

    EXPECT(array::last_extent(view) == nproma);
    EXPECT(array::conforms_layout_right(view));
    EXPECT(is_last_dimension_aligned(view, 64));

    auto span = make_mdspan<fixed_last_dim<nproma>, layout_right,
                            restrict_aligned_accessor_policy<64>>(view);

    static_assert(decltype(span)::extents_type::static_extent(2) == nproma);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type,
                                 restrict_aligned_accessor<double, 64>>);
    EXPECT(span.data_handle() == data);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.extent(1) == 4);
    EXPECT(span.extent(2) == nproma);
    EXPECT(span(1, 3, 31) == 71.);
}

CASE("test_make_mdspan_documented_simd_fast_path_example_compiles") {
    // The documented layout/alignment dispatch selects a fixed-extent SIMD kernel correctly.
    constexpr idx_t nproma = 32;
    alignas(64) double data[2 * 4 * nproma]{};
    auto view = LocalView<double, 3>{data, ArrayShape{2, 4, nproma}};

    view(1, 3, 31) = 73.;

    auto kernel = [](auto span) {
        return span(1, 3, 31);
    };

    auto result = [&]() {
        if (array::conforms_layout<layout_right>(view) && is_last_dimension_aligned(view, 64)) {
            using layout = layout_right;
            using accessor_policy = restrict_aligned_accessor_policy<64>;
            using Span32 = decltype(make_mdspan<fixed_last_dim<32>, layout, accessor_policy>(view));
            static_assert(Span32::extents_type::static_extent(2) == 32);
            static_assert(std::is_same_v<typename Span32::layout_type, layout_right>);
            static_assert(std::is_same_v<typename Span32::accessor_type, restrict_aligned_accessor<double, 64>>);
            switch (array::last_extent(view)) {
                case 16: return kernel(make_mdspan<fixed_last_dim<16>, layout, accessor_policy>(view));
                case 32: return kernel(make_mdspan<fixed_last_dim<32>, layout, accessor_policy>(view));
                case 64: return kernel(make_mdspan<fixed_last_dim<64>, layout, accessor_policy>(view));
                default: return kernel(make_mdspan<layout, accessor_policy>(view));
            }
        }
        return kernel(make_mdspan(view));
    }();

    EXPECT(result == 73.);
}

CASE("test_make_mdspan_vector_default_1d") {
    // Default vector conversion creates a rank-1 row-major span over the original storage.
    std::vector<double> data{1., 2., 3., 4.};

    auto span = make_mdspan(data);

    static_assert(decltype(span)::rank() == 1);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, default_accessor<double>>);
    EXPECT(span.extent(0) == data.size());
    EXPECT(span.data_handle() == data.data());
    EXPECT(span(2) == 3.);
}

CASE("test_make_mdspan_vector_fixed_inner_extent_const") {
    // Constness changes element mutability, not the rank or fixed inner extent of vector-of-array storage.
    std::vector<std::array<double, 3>> values{{1., 2., 3.}, {4., 5., 6.}};
    auto writable = make_mdspan(values);
    auto readonly = make_mdspan(std::as_const(values));

    static_assert(std::is_same_v<typename decltype(writable)::accessor_type, default_accessor<double>>);
    static_assert(std::is_same_v<typename decltype(readonly)::accessor_type, default_accessor<const double>>);
    static_assert(std::is_same_v<decltype(writable), decltype(make_mdspan<layout_right>(values))>);
    static_assert(std::is_same_v<decltype(readonly), decltype(make_mdspan<layout_right>(std::as_const(values)))>);

    auto check_policies = [](auto& container) {
        using Element = std::conditional_t<std::is_const_v<std::remove_reference_t<decltype(container)>>, const double, double>;
        auto check = [&](auto span) {
            using Span = decltype(span);
            static_assert(Span::rank() == 2);
            static_assert(Span::extents_type::static_extent(1) == 3);
            static_assert(std::is_same_v<typename Span::element_type, Element>);
            EXPECT(span.data_handle() == container.front().data());
            EXPECT(span.extent(0) == 2);
            EXPECT(span.extent(1) == 3);
            EXPECT(span.stride(0) == 3);
            EXPECT(span.stride(1) == 1);
            EXPECT(span(1, 2) == 6.);
        };
        check(make_mdspan<layout_right>(container));
        check(make_mdspan<layout_stride>(container));
        check(make_mdspan<default_accessor>(container));
        check(make_mdspan<aligned_accessor>(container));
        check(make_mdspan<layout_stride, default_accessor>(container));
        check(make_mdspan<default_accessor<Element>>(container));
        check(make_mdspan<layout_stride, default_accessor<Element>>(container));
        check(make_mdspan<aligned_accessor_policy<alignof(double)>>(container));
    };
    check_policies(values);
    check_policies(std::as_const(values));

    static_assert(decltype(writable)::rank() == 2);
    static_assert(decltype(readonly)::rank() == 2);
    static_assert(decltype(writable)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(readonly)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(writable)::extents_type::static_extent(1) == 3);
    static_assert(decltype(readonly)::extents_type::static_extent(1) == 3);
    static_assert(std::is_same_v<typename decltype(writable)::element_type, double>);
    static_assert(std::is_same_v<typename decltype(readonly)::element_type, const double>);
    static_assert(std::is_assignable_v<decltype(writable(1, 2)), double>);
    static_assert(!std::is_assignable_v<decltype(readonly(1, 2)), double>);

    EXPECT(writable.data_handle() == values.front().data());
    EXPECT(readonly.data_handle() == values.front().data());
    EXPECT(writable.extent(0) == values.size());
    EXPECT(readonly.extent(0) == values.size());
    EXPECT(writable.extent(1) == 3);
    EXPECT(readonly.extent(1) == 3);
    EXPECT(writable.stride(0) == 3);
    EXPECT(readonly.stride(0) == 3);
    EXPECT(writable.stride(1) == 1);
    EXPECT(readonly.stride(1) == 1);
    EXPECT(readonly(0, 1) == 2.);
    EXPECT(readonly(1, 2) == 6.);
    writable(1, 2) = 42.;
    EXPECT(values[1][2] == 42.);
    EXPECT(readonly(1, 2) == 42.);
}

CASE("test_make_mdspan_empty_contiguous_strides") {
    std::vector<double> values;
    auto check = [](const auto& span, const auto& shape, const auto& expected_strides) {
        EXPECT(span.size() == 0);
        EXPECT(span.mapping().required_span_size() == 0);
        for (std::size_t dimension = 0; dimension < span.rank(); ++dimension) {
            EXPECT(span.extent(dimension) == shape[dimension]);
            EXPECT(span.stride(dimension) == expected_strides[dimension]);
            EXPECT(span.stride(dimension) > 0);
        }
    };
    auto check_container = [&](auto& container) {
        for (const auto& shape : {std::array<std::size_t, 2>{2, 0},
                                  std::array<std::size_t, 2>{0, 3},
                                  std::array<std::size_t, 2>{0, 0}}) {
            const std::array<std::size_t, 2> strides{shape[1] == 0 ? 1 : shape[1], 1};
            auto span = make_mdspan<layout_stride>(container, shape);
            EXPECT(span.data_handle() == container.data());
            check(span, shape, strides);
        }
        std::array<std::size_t, 3> shape{2, 0, 3};
        check(make_mdspan<layout_stride>(container, shape), shape, std::array<std::size_t, 3>{3, 3, 1});
    };
    check_container(values);
    check_container(std::as_const(values));
    std::array<std::size_t, 2> shape{2, 0};
    check(make_mdspan<layout_stride>(values.data(), shape), shape, std::array<std::size_t, 2>{1, 1});
}

CASE("test_make_mdspan_explicit_extents_matching_rank") {
    // Mixed static/dynamic extents retain the source rank and shape in both conversion forms.
    double data[6]{1., 2., 3., 4., 5., 6.};
    auto view = make_local_view(data);
    using Extents = extents<std::size_t, 2, dynamic_extent>;
    auto span = make_mdspan<Extents, layout_right>(view);
    auto shaped = make_mdspan<Extents, layout_right>(view, std::array<std::size_t, 2>{2, 3});
    static_assert(decltype(span)::rank() == 2);
    static_assert(decltype(shaped)::rank() == 2);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.extent(1) == 3);
    EXPECT(span(1, 2) == 6.);
    EXPECT(shaped.extent(0) == 2);
    EXPECT(shaped.extent(1) == 3);
    EXPECT(shaped(1, 2) == 6.);
}

CASE("test_make_mdspan_vector_array_reshape") {
    // Shaped vector-of-array conversions use scalar storage with every accessor spelling.
    std::vector<std::array<double, 3>> values{{1., 2., 3.}, {4., 5., 6.}};
    auto check_policies = [](auto& container) {
        using Element = std::conditional_t<std::is_const_v<std::remove_reference_t<decltype(container)>>, const double, double>;
        auto check = [&](auto span) {
            using Span = decltype(span);
            static_assert(Span::rank() == 2);
            static_assert(std::is_same_v<typename Span::element_type, Element>);
            static_assert(std::is_assignable_v<decltype(span(1, 2)), double> == !std::is_const_v<Element>);
            EXPECT(span.data_handle() == container.front().data());
            EXPECT(span.extent(0) == 2);
            EXPECT(span.extent(1) == 3);
            EXPECT(span.stride(0) == 3);
            EXPECT(span.stride(1) == 1);
            EXPECT(span(1, 2) == 6.);
        };
        std::array<std::size_t, 2> shape{2, 3};
        check(make_mdspan(container, shape));
        check(make_mdspan<layout_stride>(container, shape));
        check(make_mdspan<default_accessor>(container, shape));
        check(make_mdspan<aligned_accessor>(container, shape));
        check(make_mdspan<layout_stride, default_accessor>(container, shape));
        check(make_mdspan<default_accessor<Element>>(container, shape));
        check(make_mdspan<layout_stride, default_accessor<Element>>(container, shape));
        check(make_mdspan<aligned_accessor_policy<alignof(double)>>(container, shape));
        check(make_mdspan<layout_stride, aligned_accessor_policy<alignof(double)>>(container, shape));
        extents<std::size_t, dynamic_extent, 3> fixed_shape{2};
        auto fixed = make_mdspan(container, fixed_shape);
        static_assert(decltype(fixed)::extents_type::static_extent(1) == 3);
        check(fixed);
    };
    check_policies(values);
    check_policies(std::as_const(values));
    auto writable = make_mdspan(values, std::array<std::size_t, 2>{2, 3});
    auto readonly = make_mdspan(std::as_const(values), std::array<std::size_t, 2>{2, 3});
    writable(1, 2) = 42.;
    EXPECT(values[1][2] == 42.);
    EXPECT(readonly(1, 2) == 42.);
}

CASE("test_make_mdspan_vector_reshape") {
    // A runtime shape exposes flat vector storage as a rank-2 row-major span.
    std::vector<double> data{1., 2., 3., 4., 5., 6.};

    std::array<std::size_t, 2> shape{2, 3};
    auto span = make_mdspan(data, shape);

    static_assert(decltype(span)::rank() == 2);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.extent(1) == 3);
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_vector_reshape_extents_layout_accessor") {
    // Vector reshaping honors mixed dynamic/fixed extents, explicit layout, and accessor selection.
    std::vector<double> data{1., 2., 3., 4., 5., 6.};

    extents<std::size_t, dynamic_extent, 3> shape{2};
    auto span = make_mdspan<layout_stride, default_accessor>(data, shape);

    static_assert(decltype(span)::rank() == 2);
    static_assert(decltype(span)::extents_type::static_extent(1) == 3);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_stride>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, default_accessor<double>>);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.extent(1) == 3);
    EXPECT(span.stride(0) == 3);
    EXPECT(span.stride(1) == 1);
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_pointer_reshape") {
    // Raw storage plus a runtime shape produces a non-owning rank-2 span with default policies.
    double data[6]{1., 2., 3., 4., 5., 6.};

    std::array<std::size_t, 2> shape{2, 3};
    auto span = make_mdspan(data, shape);

    static_assert(decltype(span)::rank() == 2);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, default_accessor<double>>);
    EXPECT(span.data_handle() == data);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.extent(1) == 3);
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_pointer_reshape_explicit_dims") {
    // Explicit dims<2> permits a braced runtime shape and preserves the requested extents type.
    double data[6]{1., 2., 3., 4., 5., 6.};

    auto span = make_mdspan<dims<2>>(data, {2, 3});

    static_assert(decltype(span)::rank() == 2);
    static_assert(std::is_same_v<typename decltype(span)::extents_type, dims<2>>);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, default_accessor<double>>);
    EXPECT(span.data_handle() == data);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.extent(1) == 3);
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_pointer_reshape_extents_layout_accessor") {
    // Raw-pointer conversion supports fixed inner extents with explicit stride layout and accessor.
    double data[6]{1., 2., 3., 4., 5., 6.};

    extents<std::size_t, dynamic_extent, 3> shape{2};
    auto span = make_mdspan<layout_stride, default_accessor>(data, shape);

    static_assert(decltype(span)::rank() == 2);
    static_assert(decltype(span)::extents_type::static_extent(1) == 3);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_stride>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, default_accessor<double>>);
    EXPECT(span.data_handle() == data);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.extent(1) == 3);
    EXPECT(span.stride(0) == 3);
    EXPECT(span.stride(1) == 1);
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_pointer_reshape_accessor_first") {
    // Accessor-first policy syntax selects an aligned accessor for reshaped raw storage.
    double data[6]{1., 2., 3., 4., 5., 6.};

    std::array<std::size_t, 2> shape{2, 3};
    auto span = make_mdspan<aligned_accessor_policy<alignof(double)>>(data, shape);

    static_assert(decltype(span)::rank() == 2);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, aligned_accessor<double, alignof(double)>>);
    EXPECT(span.data_handle() == data);
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_dynamic_fixed_extents") {
    // A 4x3 ArrayView retains its dynamic outer size when the inner extent is fixed in the type.
    std::unique_ptr<Array> array{Array::create<double>(4, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(3, 1) = 21.;

    auto span = make_mdspan<DynamicFirstFixedLastExtents, layout_stride>(view);

    expect_dynamic_first_fixed_last_extents<decltype(span)>();
    expect_same_data(span, view);
    expect_shape_2d(span, view);
    EXPECT(span(3, 1) == 21.);
}

CASE("test_make_mdspan_dynamic_fixed_extents_argument") {
    std::unique_ptr<Array> array{Array::create<double>(4, 3)};
    // A shaped 4x3 conversion uses the supplied dynamic extent and preserves the fixed inner extent.
    auto view = make_host_view<double, 2>(*array);

    view(3, 1) = 24.;

    DynamicFirstFixedLastExtents extents{4};
    auto span = make_mdspan<layout_stride>(view, extents);

    expect_dynamic_first_fixed_last_extents<decltype(span)>();
    expect_same_data(span, view);
    expect_shape_2d(span, view);
    EXPECT(span(3, 1) == 24.);
}

CASE("test_make_mdspan_localview") {
    // LocalView conversion preserves the original data, shape, and stride mapping.
    double data[6]{};
    auto view = make_local_view(data);

    view(1, 2) = 9.;

    auto span = make_mdspan<layout_stride>(view);

    expect_layout_stride_2d<decltype(span)>(span, view);
    EXPECT(span(1, 2) == 9.);
}

CASE("test_make_mdspan_localview_layout_right") {
    // A contiguous LocalView supports explicit row-major conversion without a data copy.
    double data[6]{};
    auto view = make_local_view(data);

    view(1, 2) = 11.;

    auto span = make_mdspan<layout_right>(view);

    expect_layout_right_2d<decltype(span)>(span, view);
    EXPECT(span(1, 2) == 11.);
}

CASE("test_make_mdspan_const_localview") {
    // A const LocalView object retains mutable element type and the same storage mapping.
    double data[6]{};
    auto view = make_local_view(data);
    view(0, 2) = 13.;

    const auto& const_view = view;
    auto span = make_mdspan<layout_stride>(const_view);

    static_assert(std::is_same_v<typename decltype(span)::element_type, double>);
    expect_layout_stride_2d<decltype(span)>(span, const_view);
    EXPECT(span(0, 2) == 13.);
}

CASE("test_make_mdspan_localview_runtime_extents") {
    // Supplying runtime extents to LocalView conversion preserves its shape and strides.
    double data[6]{};
    auto view = make_local_view(data);

    std::array<std::size_t, 2> extents{2, 3};
    auto span = make_mdspan<layout_stride>(view, extents);

    expect_layout_stride_2d<decltype(span)>(span, view);
}

CASE("test_make_mdspan_localview_dynamic_fixed_extents") {
    double data[6]{};
    // A 2x3 LocalView supports mixed dynamic/fixed extents with either stride or row-major layout.
    auto view = make_local_view(data);

    view(1, 2) = 17.;

    auto stride_span = make_mdspan<DynamicFirstFixedLastExtents, layout_stride>(view);
    auto right_span = make_mdspan<DynamicFirstFixedLastExtents, layout_right>(view);

    expect_dynamic_first_fixed_last_extents<decltype(stride_span)>();
    expect_layout_stride_2d<decltype(stride_span)>(stride_span, view);
    EXPECT(stride_span(1, 2) == 17.);
    expect_layout_right_2d<decltype(right_span)>(right_span, view);
    EXPECT(right_span(1, 2) == 17.);
}

CASE("test_make_mdspan_localview_dynamic_fixed_extents_argument") {
    double data[6]{};
    // A supplied extents object gives both layouts a dynamic outer size of two and fixed inner size of three.
    auto view = make_local_view(data);

    view(1, 2) = 19.;

    DynamicFirstFixedLastExtents extents{2};
    auto stride_span = make_mdspan<layout_stride>(view, extents);
    auto right_span = make_mdspan<layout_right>(view, extents);

    expect_dynamic_first_fixed_last_extents<decltype(stride_span)>();
    expect_layout_stride_2d<decltype(stride_span)>(stride_span, view);
    EXPECT(stride_span(1, 2) == 19.);
    expect_layout_right_2d<decltype(right_span)>(right_span, view);
    EXPECT(right_span(1, 2) == 19.);
}

CASE("test_make_mdspan_localview_free_function") {
    // The free conversion function exposes LocalView storage and its existing element values.
    double data[6]{};
    auto view = make_local_view(data);

    view(1, 0) = 3.;

    auto span = make_mdspan<layout_stride>(view);

    expect_same_data(span, view);
    EXPECT(span(1, 0) == 3.);
}

CASE("test_make_mdspan_localview_dynamic_fixed_extents") {
    double data[12]{};
    // A 4x3 LocalView preserves the larger runtime outer extent with a fixed inner extent.
    auto view = make_local_view_4x3(data);

    view(3, 1) = 23.;

    auto span = make_mdspan<DynamicFirstFixedLastExtents, layout_stride>(view);

    expect_dynamic_first_fixed_last_extents<decltype(span)>();
    expect_same_data(span, view);
    expect_shape_2d(span, view);
    EXPECT(span(3, 1) == 23.);
}

CASE("test_make_mdspan_localview_dynamic_fixed_extents_argument") {
    double data[12]{};
    // Explicit mixed extents preserve a 4x3 LocalView's storage and runtime shape.
    auto view = make_local_view_4x3(data);

    view(3, 1) = 29.;

    DynamicFirstFixedLastExtents extents{4};
    auto span = make_mdspan<layout_stride>(view, extents);

    expect_dynamic_first_fixed_last_extents<decltype(span)>();
    expect_same_data(span, view);
    expect_shape_2d(span, view);
    EXPECT(span(3, 1) == 29.);
}

CASE("test_make_mdspan_native_indexview") {
    // IndexView conversion preserves its index-adjusting accessor, handle, extent, and stride.
    int data[6]{};
    idx_t shape[1]{6};
    idx_t strides[1]{1};
    IndexView<int, 1> view{data, shape, strides};

    auto span = make_mdspan(view);

    static_assert(std::is_same_v<typename decltype(span)::accessor_type, index_accessor<int, ATLAS_HAVE_FORTRAN>>);
    EXPECT(span.data_handle() == view.data());
    EXPECT(span.extent(0) == 6);
    EXPECT(span.stride(0) == view.stride(0));
}

CASE("test_make_mdspan_queries_layout_right") {
    // Layout queries accept row-major views and mdspans but reject padded row strides.
    double contiguous_data[6]{};
    auto contiguous_view = make_local_view(contiguous_data);

    double strided_data[8]{};
    idx_t shape[2]{2, 3};
    idx_t padded_strides[2]{4, 1};
    auto padded_view = LocalView<double, 2>{strided_data, shape, padded_strides};

    int index_data[6]{};
    idx_t index_shape[1]{6};
    idx_t index_strides[1]{1};
    auto index_view = IndexView<int, 1>{index_data, index_shape, index_strides};

    auto stride_span = make_mdspan<layout_stride>(contiguous_view);
    auto right_span = make_mdspan<layout_right>(contiguous_view);

    EXPECT(conforms_layout_right(LayoutRightViewWithoutStrides{}));
    EXPECT(conforms_layout_right(contiguous_view));
    EXPECT(conforms_layout_right(stride_span));
    EXPECT(conforms_layout_right(right_span));
    EXPECT(conforms_layout_right(index_view));
    EXPECT(not conforms_layout_right(padded_view));
    EXPECT(array::last_extent(contiguous_view) == 3);
    EXPECT(array::last_extent(stride_span) == 3);
    EXPECT(array::last_extent(contiguous_view) != 4);
}

CASE("test_make_mdspan_queries_alignment") {
    // Base, stride, and slice-start alignment are distinct; zero alignment is rejected.
    alignas(64) double aligned_data[12]{};
    idx_t shape[2]{2, 3};
    idx_t contiguous_strides[2]{3, 1};
    idx_t padded_strides[2]{4, 1};
    auto contiguous_view = LocalView<double, 2>{aligned_data, shape, contiguous_strides};
    auto padded_view = LocalView<double, 2>{aligned_data, shape, padded_strides};
    auto span = make_mdspan<layout_stride>(padded_view);

    EXPECT(is_aligned(contiguous_view, 64));
    EXPECT(is_aligned(span, 64));
    EXPECT(not is_aligned(contiguous_view, 0));
    EXPECT(is_dimension_aligned<0>(contiguous_view, 64));
    EXPECT(not is_stride_aligned<0>(contiguous_view, 64));
    EXPECT(is_stride_aligned<0>(padded_view, 32));
    EXPECT(not is_stride_aligned<1>(padded_view, 32));
    EXPECT(not is_stride_aligned<0>(padded_view, 0));
    EXPECT(not is_dimension_aligned<0>(padded_view, 0));
    EXPECT(is_dimension_aligned<0>(contiguous_view, alignof(double)));
    EXPECT(not is_dimension_aligned<1>(contiguous_view, 64));
    EXPECT(not is_last_dimension_aligned(contiguous_view, 64));
    EXPECT(is_last_dimension_aligned(padded_view, 32));
    EXPECT(is_last_dimension_aligned(span, 32));
}

CASE("test_make_mdspan_alignment_all_slice_starts") {
    // Slice alignment checks every varying outer stride, skips singletons, and supports const inputs.
    alignas(64) double data[256]{};
    idx_t shape[3]{2, 2, 3};
    idx_t strides[3]{33, 8, 1};
    auto view = LocalView<double, 3>{data, shape, strides};
    auto span = make_mdspan<layout_stride>(view);

    EXPECT(is_aligned(view, 64));
    EXPECT(is_stride_aligned<1>(view, 64));
    EXPECT(not is_stride_aligned<0>(view, 64));
    EXPECT(is_dimension_aligned<0>(view, 64));
    EXPECT(not is_dimension_aligned<1>(view, 64));
    EXPECT(not is_last_dimension_aligned(view, 64));
    EXPECT(not is_last_dimension_aligned(span, 64));

    idx_t aligned_strides[3]{32, 8, 1};
    auto aligned_view = LocalView<double, 3>{data, shape, aligned_strides};
    EXPECT(is_last_dimension_aligned(aligned_view, 64));
    const auto& const_aligned_view = aligned_view;
    EXPECT(is_aligned(const_aligned_view, 64));
    EXPECT(is_aligned<64>(const_aligned_view));
    EXPECT(is_dimension_aligned<2>(const_aligned_view, 64));
    EXPECT(is_last_dimension_aligned(const_aligned_view, 64));
    EXPECT(is_aligned(make_mdspan(aligned_view), 64));
    EXPECT(is_aligned<64>(make_mdspan(aligned_view)));
    EXPECT(is_stride_aligned<0>(make_mdspan(aligned_view), 64));
    EXPECT(is_dimension_aligned<2>(make_mdspan(aligned_view), 64));
    EXPECT(is_last_dimension_aligned(make_mdspan(aligned_view), 64));
    for (idx_t outer = 0; outer < shape[0]; ++outer) {
        for (idx_t middle = 0; middle < shape[1]; ++middle) {
            EXPECT(reinterpret_cast<std::uintptr_t>(&aligned_view(outer, middle, 0)) % 64 == 0);
        }
    }

    idx_t singleton_shape[3]{1, 2, 3};
    auto singleton_view = LocalView<double, 3>{data, singleton_shape, strides};
    EXPECT(not is_stride_aligned<0>(singleton_view, 64));
    EXPECT(is_last_dimension_aligned(singleton_view, 64));

    idx_t vector_shape[1]{3};
    idx_t vector_strides[1]{1};
    auto vector_view = LocalView<double, 1>{data, vector_shape, vector_strides};
    EXPECT(is_last_dimension_aligned(vector_view, 64));
    EXPECT(not is_stride_aligned<0>(vector_view, 64));

    auto offset_view = LocalView<double, 3>{data + 1, shape, aligned_strides};
    const auto& const_offset_view = offset_view;
    EXPECT(is_stride_aligned<0>(const_offset_view, 64));
    EXPECT(not is_aligned(const_offset_view, 64));
    EXPECT(not is_dimension_aligned<0>(const_offset_view, 64));
    EXPECT(not is_last_dimension_aligned(const_offset_view, 64));
}

CASE("test_make_mdspan_queries_unaligned") {
    // An offset base fails stronger alignment requirements even when its strides are contiguous.
    alignas(16) double data[7]{};
    idx_t shape[2]{2, 3};
    idx_t strides[2]{3, 1};
    auto view = LocalView<double, 2>{data + 1, shape, strides};

    EXPECT(not is_aligned(view, 16));
    EXPECT(not is_dimension_aligned<1>(view, 16));
    EXPECT(not is_last_dimension_aligned(view, 16));
    EXPECT(is_aligned(view, alignof(double)));
}

CASE("test_aligned_accessor") {
    // Aligned accessors expose their alignment and allow only supported const/alignment conversions.
    alignas(64) double data[6]{};
    using Span = mdspan<double, extents<std::size_t, 2, 3>, layout_right, aligned_accessor<double, 64>>;
    Span span{data};

    span(1, 2) = 43.;

    static_assert(std::is_same_v<typename Span::accessor_type, aligned_accessor<double, 64>>);
    static_assert(aligned_accessor<double, 64>::byte_alignment == 64);
    static_assert(std::is_convertible_v<aligned_accessor<double, 64>, aligned_accessor<const double, 32>>);
    static_assert(std::is_convertible_v<aligned_accessor<double, 64>, default_accessor<const double>>);
    static_assert(not std::is_convertible_v<default_accessor<double>, aligned_accessor<double, 64>>);
    static_assert(std::is_constructible_v<aligned_accessor<double, 64>, default_accessor<double>>);

    EXPECT(span.data_handle() == data);
    EXPECT(span(1, 2) == 43.);
    EXPECT(is_sufficiently_aligned<64>(data));
    EXPECT(not is_sufficiently_aligned<64>(data + 1));
}

CASE("test_restrict_aligned_accessor") {
    // Restrict-aligned accessors retain aligned access and support compatible accessor conversions.
    alignas(64) double data[6]{};
    using Span = mdspan<double, extents<std::size_t, 2, 3>, layout_right, restrict_aligned_accessor<double, 64>>;
    Span span{data};

    span(0, 2) = 47.;

    static_assert(std::is_same_v<typename Span::accessor_type, restrict_aligned_accessor<double, 64>>);
    static_assert(restrict_aligned_accessor<double, 64>::byte_alignment == 64);
    static_assert(std::is_convertible_v<restrict_aligned_accessor<double, 64>, restrict_aligned_accessor<const double, 32>>);
    static_assert(std::is_convertible_v<restrict_aligned_accessor<double, 64>, restrict_accessor<const double>>);
    static_assert(std::is_convertible_v<restrict_aligned_accessor<double, 64>, aligned_accessor<const double, 64>>);
    static_assert(std::is_convertible_v<restrict_aligned_accessor<double, 64>, default_accessor<const double>>);
    static_assert(std::is_constructible_v<restrict_aligned_accessor<double, 64>, aligned_accessor<double, 64>>);
    static_assert(std::is_constructible_v<restrict_aligned_accessor<double, 64>, restrict_accessor<double>>);
    static_assert(std::is_constructible_v<restrict_aligned_accessor<double, 64>, default_accessor<double>>);

    EXPECT(span.data_handle() == data);
    EXPECT(span(0, 2) == 47.);
}

CASE("test_make_mdspan_accessor_first_syntax") {
    // Accessor-first syntax selects each accessor family while preserving writable shared storage.
    double data[6]{};
    auto view = make_local_view(data);

    auto default_span = make_mdspan<default_accessor>(view);
    auto aligned_span = make_mdspan<aligned_accessor_policy<alignof(double)>>(view);
    auto restrict_span = make_mdspan<restrict_accessor>(view);
    auto restrict_aligned_span = make_mdspan<restrict_aligned_accessor_policy<alignof(double)>>(view);

    static_assert(std::is_same_v<typename decltype(default_span)::accessor_type, default_accessor<double>>);
    static_assert(std::is_same_v<typename decltype(aligned_span)::accessor_type, aligned_accessor<double, alignof(double)>>);
    static_assert(std::is_same_v<typename decltype(restrict_span)::accessor_type, restrict_accessor<double>>);
    static_assert(std::is_same_v<typename decltype(restrict_aligned_span)::accessor_type, restrict_aligned_accessor<double, alignof(double)>>);

    default_span(1, 0) = 51.;
    aligned_span(1, 1) = 52.;
    restrict_span(1, 2) = 53.;
    restrict_aligned_span(0, 2) = 54.;

    EXPECT(default_span.data_handle() == view.data());
    EXPECT(aligned_span.data_handle() == view.data());
    EXPECT(restrict_span.data_handle() == view.data());
    EXPECT(restrict_aligned_span.data_handle() == view.data());
    EXPECT(view(1, 0) == 51.);
    EXPECT(view(1, 1) == 52.);
    EXPECT(view(1, 2) == 53.);
    EXPECT(view(0, 2) == 54.);
}

CASE("test_make_mdspan_accessor_first_explicit_alignment") {
    // Explicit alignment policies work with both inferred and supplied shapes.
    alignas(64) double data[6]{};
    auto view = make_local_view(data);

    std::array<std::size_t, 2> shape{2, 3};
    auto aligned_span = make_mdspan<aligned_accessor_policy<64>>(view);
    auto restrict_aligned_span = make_mdspan<restrict_aligned_accessor_policy<64>>(view, shape);

    static_assert(std::is_same_v<typename decltype(aligned_span)::accessor_type, aligned_accessor<double, 64>>);
    static_assert(std::is_same_v<typename decltype(restrict_aligned_span)::accessor_type, restrict_aligned_accessor<double, 64>>);
    static_assert(std::is_same_v<typename decltype(aligned_span)::layout_type, layout_stride>);
    static_assert(std::is_same_v<typename decltype(restrict_aligned_span)::layout_type, layout_stride>);

    aligned_span(1, 1) = 65.;
    restrict_aligned_span(1, 2) = 66.;

    EXPECT(aligned_span.data_handle() == view.data());
    EXPECT(restrict_aligned_span.data_handle() == view.data());
    EXPECT(aligned_span(1, 1) == 65.);
    EXPECT(restrict_aligned_span(1, 2) == 66.);
    EXPECT(is_sufficiently_aligned<64>(data));
}

CASE("test_make_mdspan_shaped_accessor_uses_requested_extents") {
    // Accessor-first shaping must replace source static extents with the requested output extents.
    double data[6]{1., 2., 3., 4., 5., 6.};
    using Source = mdspan<double, extents<std::size_t, 6>, layout_right>;
    Source source{data};
    std::array<std::size_t, 1> shape{3};
    auto ordinary = make_mdspan(source, shape);
    auto aligned = make_mdspan<aligned_accessor>(source, shape);
    auto restrict_aligned = make_mdspan<restrict_aligned_accessor>(source, shape);
    auto temporary = make_mdspan<aligned_accessor>(Source{data}, shape);
    auto fixed = make_mdspan<aligned_accessor>(source, extents<std::size_t, 3>{});

    static_assert(decltype(aligned)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(restrict_aligned)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(fixed)::extents_type::static_extent(0) == 3);
    static_assert(std::is_same_v<typename decltype(aligned)::accessor_type,
                                 aligned_accessor<double, alignof(double)>>);
    EXPECT(source.extent(0) == 6);
    EXPECT(ordinary.extent(0) == 3);
    EXPECT(aligned.extent(0) == ordinary.extent(0));
    EXPECT(restrict_aligned.extent(0) == 3);
    EXPECT(temporary.extent(0) == 3);
    EXPECT(fixed.extent(0) == 3);
    EXPECT(aligned.data_handle() == data);
    EXPECT(aligned.stride(0) == 1);
    EXPECT(aligned(2) == ordinary(2));
    EXPECT(restrict_aligned(2) == 3.);
    EXPECT(temporary(2) == 3.);
    EXPECT(fixed(2) == 3.);
}

CASE("test_make_mdspan_accessor_type") {
    // Bound accessor types work for views and pointers, including explicitly read-only elements.
    alignas(64) double data[6]{};
    auto view = make_local_view(data);

    std::array<std::size_t, 2> shape{2, 3};
    auto aligned_span = make_mdspan<aligned_accessor<double, 64>>(view);
    auto restrict_aligned_span = make_mdspan<layout_stride, restrict_aligned_accessor<double, 64>>(view, shape);
    auto pointer_span = make_mdspan<aligned_accessor<double, 64>>(data, shape);
    auto const_pointer_span = make_mdspan<default_accessor<const double>>(data, shape);

    static_assert(std::is_same_v<typename decltype(aligned_span)::accessor_type, aligned_accessor<double, 64>>);
    static_assert(std::is_same_v<typename decltype(restrict_aligned_span)::accessor_type, restrict_aligned_accessor<double, 64>>);
    static_assert(std::is_same_v<typename decltype(pointer_span)::accessor_type, aligned_accessor<double, 64>>);
    static_assert(std::is_same_v<typename decltype(const_pointer_span)::accessor_type, default_accessor<const double>>);
    static_assert(std::is_const_v<typename decltype(const_pointer_span)::element_type>);

    aligned_span(1, 0) = 67.;
    restrict_aligned_span(1, 1) = 68.;
    pointer_span(1, 2) = 69.;

    EXPECT(aligned_span.data_handle() == view.data());
    EXPECT(restrict_aligned_span.data_handle() == view.data());
    EXPECT(pointer_span.data_handle() == data);
    EXPECT(const_pointer_span.data_handle() == data);
    EXPECT(view(1, 0) == 67.);
    EXPECT(view(1, 1) == 68.);
    EXPECT(view(1, 2) == 69.);
}

CASE("test_make_mdspan_container_accessor_type") {
    // A reshaped container accepts a bound accessor type and explicit row-major layout.
    std::vector<double> values{1., 2., 3., 4., 5., 6.};

    std::array<std::size_t, 2> shape{2, 3};
    auto span = make_mdspan<layout_right, default_accessor<double>>(values, shape);

    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, default_accessor<double>>);
    EXPECT(span.data_handle() == values.data());
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_from_mdspan_changes_layout") {
    // Converting an existing mdspan can change layout while preserving accessor, shape, and storage.
    double data[6]{};
    auto view = make_local_view(data);
    auto stride_span = make_mdspan<layout_stride>(view);

    stride_span(1, 2) = 31.;
    auto right_span = make_mdspan<layout_right>(stride_span);

    static_assert(std::is_same_v<typename decltype(right_span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(right_span)::accessor_type, default_accessor<double>>);
    EXPECT(right_span.data_handle() == stride_span.data_handle());
    EXPECT(right_span.extent(0) == stride_span.extent(0));
    EXPECT(right_span.extent(1) == stride_span.extent(1));
    EXPECT(right_span(1, 2) == 31.);
}

CASE("test_make_mdspan_from_mdspan_changes_accessor") {
    // Converting an existing mdspan can replace its accessor and layout without copying elements.
    double data[6]{};
    auto view = make_local_view(data);
    auto restrict_span = make_mdspan<layout_right, restrict_accessor>(view);

    auto preserved = make_mdspan(restrict_span);
    auto preserved_shaped = make_mdspan(restrict_span, std::array<std::size_t, 2>{2, 3});
    static_assert(std::is_same_v<typename decltype(preserved)::accessor_type, restrict_accessor<double>>);
    static_assert(std::is_same_v<typename decltype(preserved_shaped)::accessor_type, restrict_accessor<double>>);

    restrict_span(0, 2) = 37.;
    EXPECT(preserved.data_handle() == restrict_span.data_handle());
    EXPECT(preserved_shaped.data_handle() == restrict_span.data_handle());
    EXPECT(preserved(0, 2) == 37.);
    EXPECT(preserved_shaped(0, 2) == 37.);
    auto default_span = make_mdspan<layout_stride, default_accessor>(restrict_span);

    static_assert(std::is_same_v<typename decltype(default_span)::layout_type, layout_stride>);
    static_assert(std::is_same_v<typename decltype(default_span)::accessor_type, default_accessor<double>>);
    EXPECT(default_span.data_handle() == restrict_span.data_handle());
    EXPECT(default_span.extent(0) == restrict_span.extent(0));
    EXPECT(default_span.extent(1) == restrict_span.extent(1));
    EXPECT(default_span(0, 2) == 37.);
}

CASE("test_make_mdspan_from_mdspan_fixed_extents") {
    // Existing mdspan conversion can refine the inner extent to a compile-time constant.
    double data[6]{};
    auto view = make_local_view(data);
    auto stride_span = make_mdspan<layout_stride>(view);

    stride_span(1, 2) = 41.;
    auto fixed_span = make_mdspan<DynamicFirstFixedLastExtents, layout_right>(stride_span);

    expect_dynamic_first_fixed_last_extents<decltype(fixed_span)>();
    static_assert(std::is_same_v<typename decltype(fixed_span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(fixed_span)::accessor_type, default_accessor<double>>);
    EXPECT(fixed_span.data_handle() == stride_span.data_handle());
    EXPECT(fixed_span.extent(0) == 2);
    EXPECT(fixed_span.extent(1) == 3);
    EXPECT(fixed_span(1, 2) == 41.);
}

CASE("test_make_mdspan_from_mdspan_preserves_extents") {
    // Fully static extents survive default, layout-changing, accessor-changing, and shaped conversions.
    double data[6]{};
    using FixedExtents = extents<std::size_t, 2, 3>;
    using InputSpan = mdspan<double, FixedExtents, layout_right, aligned_accessor<double, alignof(double)>>;
    InputSpan fixed_span{data};

    fixed_span(1, 2) = 61.;
    auto same_span = make_mdspan(fixed_span);
    auto stride_span = make_mdspan<layout_stride>(fixed_span);
    auto default_span = make_mdspan<default_accessor>(fixed_span);
    auto shaped_span = make_mdspan(fixed_span, FixedExtents{});

    static_assert(std::is_same_v<typename decltype(same_span)::extents_type, FixedExtents>);
    static_assert(std::is_same_v<typename decltype(stride_span)::extents_type, FixedExtents>);
    static_assert(std::is_same_v<typename decltype(default_span)::extents_type, FixedExtents>);
    static_assert(std::is_same_v<typename decltype(shaped_span)::extents_type, FixedExtents>);
    static_assert(std::is_same_v<typename decltype(same_span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(stride_span)::layout_type, layout_stride>);
    static_assert(std::is_same_v<typename decltype(default_span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(shaped_span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(same_span)::accessor_type, typename InputSpan::accessor_type>);
    static_assert(std::is_same_v<typename decltype(stride_span)::accessor_type, typename InputSpan::accessor_type>);
    static_assert(std::is_same_v<typename decltype(default_span)::accessor_type, default_accessor<double>>);
    static_assert(std::is_same_v<typename decltype(shaped_span)::accessor_type, typename InputSpan::accessor_type>);
    EXPECT(same_span.data_handle() == fixed_span.data_handle());
    EXPECT(stride_span.data_handle() == fixed_span.data_handle());
    EXPECT(default_span.data_handle() == fixed_span.data_handle());
    EXPECT(shaped_span.data_handle() == fixed_span.data_handle());
    EXPECT(same_span(1, 2) == 61.);
    EXPECT(stride_span(1, 2) == 61.);
    EXPECT(default_span(1, 2) == 61.);
    EXPECT(shaped_span(1, 2) == 61.);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
