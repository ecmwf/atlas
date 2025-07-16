#include "pluto/mdspan.h"
#include <vector>
#include <array>
#include <iostream>
#include <string_view>

#if defined(__cpp_multidimensional_subscript)
#define MDINDEX(...) __VA_ARGS__
#else
#define MDINDEX(...) std::array{__VA_ARGS__}
#endif

#ifdef __GNUC__
#if __GNUC__ <= 8
#define CTAD_NOT_SUPPORTED 1
#endif
#endif

void print_features() {
    std::cout << "__cplusplus                           = " << __cplusplus << std::endl;
#if defined(__cpp_multidimensional_subscript)
    std::cout << "__cpp_multidimensional_subscript      = " << __cpp_multidimensional_subscript << std::endl;
#else
    std::cout << "__cpp_multidimensional_subscript      = not defined" << std::endl;
#endif
#if defined(__cpp_lib_mdspan)
    std::cout << "__cpp_lib_mdspan                      = " << __cpp_lib_mdspan << std::endl;
#else
    std::cout << "__cpp_lib_mdspan                      = not defined" << std::endl;
#endif
    std::cout << "PLUTO_HAVE_MDSPAN                     = " << PLUTO_HAVE_MDSPAN << std::endl;
}

template<typename MDSPAN>
void print(std::string_view name, MDSPAN span) {
    if constexpr(span.rank() == 1) {
        std::cout << name << " = {";
        for (typename MDSPAN::index_type i=0; i<span.extent(0); ++i) {
            std::cout << "\t" << span[i];
        }
        std::cout << "   }" << std::endl;
    }
    if constexpr(span.rank() == 2) {
        std::cout << name << " = {\n";
        for (typename MDSPAN::index_type i=0; i<span.extent(0); ++i) {
            for (typename MDSPAN::index_type j=0; j<span.extent(1); ++j) {
                std::cout << "\t" << span[MDINDEX(i,j)];
            }
            std::cout << "\n";
        }
        std::cout << "}" << std::endl;
    }
}
#define PRINT(span) print(#span,span);

template<typename MDSPAN, typename EXPECTED>
void expect_equal(MDSPAN span, EXPECTED expected) {
    if(span.size() != expected.size()) {
        throw std::out_of_range("size mismatch");
    }
    if constexpr(span.rank() == 1) {
        for (typename MDSPAN::index_type i=0; i<span.extent(0); ++i) {
            if(span[i] != expected[i]) {
                throw std::runtime_error("value mismatch");
            }
        }
    }
    if constexpr(span.rank() == 2) {
        int n=0;
        for (typename MDSPAN::index_type i=0; i<span.extent(0); ++i) {
            for (typename MDSPAN::index_type j=0; j<span.extent(1); ++j, ++n) {
                if(span[MDINDEX(i,j)] != expected[n]) {
                    throw std::runtime_error("mismatch");
                }
            }
        }
    }
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
    print_features();

    std::vector<double> container{1,2,3,4,5,6};

    std::cout<< "\n\nTest dynamic extents" << std::endl;
    {
#ifdef CTAD_NOT_SUPPORTED
        pluto::mdspan<double,pluto::dims<2>>                      view_2x3{container.data(),2,3};
#else
        pluto::mdspan                                             view_2x3{container.data(),2,3};
#endif
        pluto::mdspan<double,pluto::dims<2>,pluto::layout_right>  view_2x3_layout_right{container.data(),2,3};
        pluto::mdspan<double,pluto::dims<2>,pluto::layout_left>   view_2x3_layout_left{container.data(),2,3};
        pluto::mdspan<double,pluto::dims<1>,pluto::layout_stride> view_2x3_row1{container.data()+0,{pluto::dims<1>(3),std::array{1}}};
        pluto::mdspan<double,pluto::dims<1>,pluto::layout_stride> view_2x3_row2{container.data()+3,{pluto::dims<1>(3),std::array{1}}};
        pluto::mdspan<double,pluto::dims<1>,pluto::layout_stride> view_2x3_col1{container.data()+0,{pluto::dims<1>(2),std::array{3}}};
        pluto::mdspan<double,pluto::dims<1>,pluto::layout_stride> view_2x3_col2{container.data()+1,{pluto::dims<1>(2),std::array{3}}};
        pluto::mdspan<double,pluto::dims<1>,pluto::layout_stride> view_2x3_col3{container.data()+2,{pluto::dims<1>(2),std::array{3}}};
        PRINT(view_2x3);
        PRINT(view_2x3_layout_right);
        PRINT(view_2x3_layout_left);
        PRINT(view_2x3_row1);
        PRINT(view_2x3_row2);
        PRINT(view_2x3_col1);
        PRINT(view_2x3_col2);
        PRINT(view_2x3_col3);
        expect_equal(view_2x3, std::array{1,2,3,4,5,6});
        expect_equal(view_2x3_layout_right, std::array{1,2,3,4,5,6});
        expect_equal(view_2x3_layout_left, std::array{1,3,5,2,4,6});
        expect_equal(view_2x3_row1, std::array{1,2,3});
        expect_equal(view_2x3_row2, std::array{4,5,6});
        expect_equal(view_2x3_col1, std::array{1,4});
        expect_equal(view_2x3_col2, std::array{2,5});
        expect_equal(view_2x3_col3, std::array{3,6});
    }

    std::cout<< "\n\nTest static extents" << std::endl;
    {
        using shape_2x3     = pluto::extents<size_t,2,3>;
        using shape_2x3_row = pluto::extents<size_t,3>;
        using shape_2x3_col = pluto::extents<size_t,2>;
        pluto::mdspan<double,shape_2x3>                          view_2x3{container.data()};
        pluto::mdspan<double,shape_2x3,pluto::layout_right>      view_2x3_layout_right{container.data()};
        pluto::mdspan<double,shape_2x3,pluto::layout_left>       view_2x3_layout_left{container.data()};
        pluto::mdspan<double,shape_2x3_row,pluto::layout_stride> view_2x3_row1{container.data()+0,{shape_2x3_row{},std::array{1}}};
        pluto::mdspan<double,shape_2x3_row,pluto::layout_stride> view_2x3_row2{container.data()+3,{shape_2x3_row{},std::array{1}}};
        pluto::mdspan<double,shape_2x3_col,pluto::layout_stride> view_2x3_col1{container.data()+0,{shape_2x3_col{},std::array{3}}};
        pluto::mdspan<double,shape_2x3_col,pluto::layout_stride> view_2x3_col2{container.data()+1,{shape_2x3_col{},std::array{3}}};
        pluto::mdspan<double,shape_2x3_col,pluto::layout_stride> view_2x3_col3{container.data()+2,{shape_2x3_col{},std::array{3}}};
        PRINT(view_2x3);
        PRINT(view_2x3_layout_right);
        PRINT(view_2x3_layout_left);
        PRINT(view_2x3_row1);
        PRINT(view_2x3_row2);
        PRINT(view_2x3_col1);
        PRINT(view_2x3_col2);
        PRINT(view_2x3_col3);
        expect_equal(view_2x3, std::array{1,2,3,4,5,6});
        expect_equal(view_2x3_layout_right, std::array{1,2,3,4,5,6});
        expect_equal(view_2x3_layout_left, std::array{1,3,5,2,4,6});
        expect_equal(view_2x3_row1, std::array{1,2,3});
        expect_equal(view_2x3_row2, std::array{4,5,6});
        expect_equal(view_2x3_col1, std::array{1,4});
        expect_equal(view_2x3_col2, std::array{2,5});
        expect_equal(view_2x3_col3, std::array{3,6});
    }

    std::cout<< "\n\nTest assignment nonconst to const" << std::endl;
    {
        pluto::mdspan<double      ,pluto::dims<2>> view_2x3_nonconst{container.data(),2,3};
        pluto::mdspan<double const,pluto::dims<2>> view_2x3_const = view_2x3_nonconst;
        PRINT(view_2x3_const);
        expect_equal(view_2x3_const, std::array{1,2,3,4,5,6});
    }

    std::cout<< "\n\nTest layout_stride assignment" << std::endl;
    {
        using shape_2x3_row = pluto::extents<size_t,3>;
        using shape_2x3_col = pluto::extents<size_t,2>;
        pluto::mdspan<double,shape_2x3_row,pluto::layout_stride> view_2x3_row1{container.data()+0,{shape_2x3_row{},std::array{1}}};
        pluto::mdspan<double,shape_2x3_col,pluto::layout_stride> view_2x3_col1{container.data()+0,{shape_2x3_col{},std::array{3}}};

        // This should work, as default layout_right matches the stride=1
        pluto::mdspan<double,pluto::dims<1>> view_row = view_2x3_row1;
        PRINT(view_row);

        // This fails, as default layout_right doesn't match the stride=3
        //     pluto::mdspan<double,pluto::extents<size_t,2>> view_col = view_2x3_col1;
    }

}

