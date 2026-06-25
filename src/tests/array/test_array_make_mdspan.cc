#include <array>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>

#include "atlas/array/Array.h"
#include "atlas/array/LocalView.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/native/NativeIndexView.h"
#include "atlas/array/make_mdspan.h"
#include "atlas/mdspan.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;

namespace atlas {
namespace test {

using DynamicFirstFixedLastExtents = extents<std::size_t, dynamic_extent, 3>;

template <typename Span, typename View>
void expect_same_data(const Span& span, const View& view) {
    static_assert(std::is_same_v<typename Span::accessor_type, restrict_accessor<typename Span::element_type>>);
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

CASE("test_make_mdspan_arrayview_layout_stride") {
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 2) = 12.;

    auto span = make_mdspan<layout_stride>(view);

    expect_layout_stride_2d<decltype(span)>(span, view);
    EXPECT(span(1, 2) == 12.);
}

CASE("test_make_mdspan_arrayview_layout_right") {
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 2) = 42.;

    auto span = make_mdspan<layout_right>(view);

    expect_layout_right_2d<decltype(span)>(span, view);
    EXPECT(span(1, 2) == 42.);
}

CASE("test_make_mdspan_arrayview_free_function") {
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 1) = 7.;

    auto span = make_mdspan<layout_stride>(view);

    expect_same_data(span, view);
    EXPECT(span(1, 1) == 7.);
}

CASE("test_make_mdspan_const_arrayview") {
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);
    view(0, 2) = 5.;

    const auto& const_view = view;
    auto span = make_mdspan<layout_stride>(const_view);

    static_assert(std::is_const_v<typename decltype(span)::element_type>);
    expect_same_data(span, const_view);
    EXPECT(span(0, 2) == 5.);
}

CASE("test_make_mdspan_runtime_extents") {
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    std::array<std::size_t, 2> extents{2, 3};
    auto span = make_mdspan<layout_stride>(view, extents);

    expect_layout_stride_2d<decltype(span)>(span, view);
}

CASE("test_make_mdspan_arrayview_dynamic_fixed_extents") {
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
    std::unique_ptr<Array> array{Array::create<double>(2, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(1, 2) = 63.;

    auto dim_span = make_mdspan<StaticExtent<1, 3>>(view);
    auto last_span = make_mdspan<StaticLastExtent<3>, layout_right, default_accessor>(view);

    static_assert(decltype(dim_span)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(dim_span)::extents_type::static_extent(1) == 3);
    static_assert(std::is_same_v<typename decltype(dim_span)::layout_type, layout_stride>);
    static_assert(std::is_same_v<typename decltype(dim_span)::accessor_type, restrict_accessor<double>>);
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
    std::unique_ptr<Array> array{Array::create<double>(4, 3)};
    auto view = make_host_view<double, 2>(*array);

    view(3, 2) = 64.;

    std::array<std::size_t, 2> extents{4, 3};
    auto span = make_mdspan<StaticLastExtent<3>, layout_right, default_accessor>(view, extents);

    static_assert(decltype(span)::extents_type::static_extent(0) == dynamic_extent);
    static_assert(decltype(span)::extents_type::static_extent(1) == 3);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, default_accessor<double>>);
    EXPECT(span.extent(0) == 4);
    EXPECT(span.extent(1) == 3);
    EXPECT(span(3, 2) == 64.);
}

CASE("test_make_mdspan_vector_default_1d") {
    std::vector<double> data{1., 2., 3., 4.};

    auto span = make_mdspan(data);

    static_assert(decltype(span)::rank() == 1);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, restrict_accessor<double>>);
    EXPECT(span.extent(0) == data.size());
    EXPECT(span.data_handle() == data.data());
    EXPECT(span(2) == 3.);
}

CASE("test_make_mdspan_vector_reshape") {
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
    double data[6]{1., 2., 3., 4., 5., 6.};

    std::array<std::size_t, 2> shape{2, 3};
    auto span = make_mdspan(data, shape);

    static_assert(decltype(span)::rank() == 2);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, restrict_accessor<double>>);
    EXPECT(span.data_handle() == data);
    EXPECT(span.extent(0) == 2);
    EXPECT(span.extent(1) == 3);
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_pointer_reshape_extents_layout_accessor") {
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
    double data[6]{1., 2., 3., 4., 5., 6.};

    std::array<std::size_t, 2> shape{2, 3};
    auto span = make_mdspan<aligned_accessor>(data, shape);

    static_assert(decltype(span)::rank() == 2);
    static_assert(std::is_same_v<typename decltype(span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(span)::accessor_type, aligned_accessor<double, alignof(double)>>);
    EXPECT(span.data_handle() == data);
    EXPECT(span(1, 2) == 6.);
}

CASE("test_make_mdspan_dynamic_fixed_extents") {
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
    double data[6]{};
    auto view = make_local_view(data);

    view(1, 2) = 9.;

    auto span = make_mdspan<layout_stride>(view);

    expect_layout_stride_2d<decltype(span)>(span, view);
    EXPECT(span(1, 2) == 9.);
}

CASE("test_make_mdspan_localview_layout_right") {
    double data[6]{};
    auto view = make_local_view(data);

    view(1, 2) = 11.;

    auto span = make_mdspan<layout_right>(view);

    expect_layout_right_2d<decltype(span)>(span, view);
    EXPECT(span(1, 2) == 11.);
}

CASE("test_make_mdspan_const_localview") {
    double data[6]{};
    auto view = make_local_view(data);
    view(0, 2) = 13.;

    const auto& const_view = view;
    auto span = make_mdspan<layout_stride>(const_view);

    static_assert(std::is_const_v<typename decltype(span)::element_type>);
    expect_layout_stride_2d<decltype(span)>(span, const_view);
    EXPECT(span(0, 2) == 13.);
}

CASE("test_make_mdspan_localview_runtime_extents") {
    double data[6]{};
    auto view = make_local_view(data);

    std::array<std::size_t, 2> extents{2, 3};
    auto span = make_mdspan<layout_stride>(view, extents);

    expect_layout_stride_2d<decltype(span)>(span, view);
}

CASE("test_make_mdspan_localview_dynamic_fixed_extents") {
    double data[6]{};
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
    double data[6]{};
    auto view = make_local_view(data);

    view(1, 0) = 3.;

    auto span = make_mdspan<layout_stride>(view);

    expect_same_data(span, view);
    EXPECT(span(1, 0) == 3.);
}

CASE("test_make_mdspan_localview_dynamic_fixed_extents") {
    double data[12]{};
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

    EXPECT(can_use_layout_right(contiguous_view));
    EXPECT(can_use_layout_right(stride_span));
    EXPECT(can_use_layout_right(right_span));
    EXPECT(can_use_layout_right(index_view));
    EXPECT(not can_use_layout_right(padded_view));
}

CASE("test_make_mdspan_queries_alignment") {
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
    EXPECT(not is_dimension_aligned<0>(contiguous_view, 64));
    EXPECT(is_dimension_aligned<0>(contiguous_view, alignof(double)));
    EXPECT(not is_dimension_aligned<1>(contiguous_view, 64));
    EXPECT(not is_last_dimension_aligned(contiguous_view, 64));
    EXPECT(is_last_dimension_aligned(padded_view, 32));
    EXPECT(is_last_dimension_aligned(span, 32));
}

CASE("test_make_mdspan_queries_unaligned") {
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
    double data[6]{};
    auto view = make_local_view(data);

    auto default_span = make_mdspan<default_accessor>(view);
    auto aligned_span = make_mdspan<aligned_accessor>(view);
    auto restrict_span = make_mdspan<restrict_accessor>(view);
    auto restrict_aligned_span = make_mdspan<restrict_aligned_accessor>(view);

    static_assert(std::is_same_v<typename decltype(default_span)::accessor_type, default_accessor<double>>);
    static_assert(std::is_same_v<typename decltype(aligned_span)::accessor_type, aligned_accessor<double, alignof(double)>>);
    static_assert(std::is_same_v<typename decltype(restrict_span)::accessor_type, restrict_accessor<double>>);
    static_assert(std::is_same_v<typename decltype(restrict_aligned_span)::accessor_type,
                                 restrict_aligned_accessor<double, alignof(double)>>);

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

CASE("test_make_mdspan_from_mdspan_changes_layout") {
    double data[6]{};
    auto view = make_local_view(data);
    auto stride_span = make_mdspan<layout_stride>(view);

    stride_span(1, 2) = 31.;
    auto right_span = make_mdspan<layout_right>(stride_span);

    static_assert(std::is_same_v<typename decltype(right_span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(right_span)::accessor_type, restrict_accessor<double>>);
    EXPECT(right_span.data_handle() == stride_span.data_handle());
    EXPECT(right_span.extent(0) == stride_span.extent(0));
    EXPECT(right_span.extent(1) == stride_span.extent(1));
    EXPECT(right_span(1, 2) == 31.);
}

CASE("test_make_mdspan_from_mdspan_changes_accessor") {
    double data[6]{};
    auto view = make_local_view(data);
    auto restrict_span = make_mdspan<layout_right>(view);

    restrict_span(0, 2) = 37.;
    auto default_span = make_mdspan<layout_stride, default_accessor>(restrict_span);

    static_assert(std::is_same_v<typename decltype(default_span)::layout_type, layout_stride>);
    static_assert(std::is_same_v<typename decltype(default_span)::accessor_type, default_accessor<double>>);
    EXPECT(default_span.data_handle() == restrict_span.data_handle());
    EXPECT(default_span.extent(0) == restrict_span.extent(0));
    EXPECT(default_span.extent(1) == restrict_span.extent(1));
    EXPECT(default_span(0, 2) == 37.);
}

CASE("test_make_mdspan_from_mdspan_fixed_extents") {
    double data[6]{};
    auto view = make_local_view(data);
    auto stride_span = make_mdspan<layout_stride>(view);

    stride_span(1, 2) = 41.;
    auto fixed_span = make_mdspan<DynamicFirstFixedLastExtents, layout_right>(stride_span);

    expect_dynamic_first_fixed_last_extents<decltype(fixed_span)>();
    static_assert(std::is_same_v<typename decltype(fixed_span)::layout_type, layout_right>);
    static_assert(std::is_same_v<typename decltype(fixed_span)::accessor_type, restrict_accessor<double>>);
    EXPECT(fixed_span.data_handle() == stride_span.data_handle());
    EXPECT(fixed_span.extent(0) == 2);
    EXPECT(fixed_span.extent(1) == 3);
    EXPECT(fixed_span(1, 2) == 41.);
}

CASE("test_make_mdspan_from_mdspan_preserves_extents") {
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
