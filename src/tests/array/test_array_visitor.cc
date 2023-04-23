#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/helpers/ArrayVisitor.h"
#include "atlas/library/config.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;
using namespace atlas::array::helpers;

namespace atlas {
namespace test {

    CASE("test_array_slicer_2d") {

    ArrayT<double> arr1(3, 5);
    auto view1 = make_view<double, 2>(arr1);
    view1.assign({
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

    ArrayT<double> arr2(3, 5, 1);
    auto view2 = make_view<double, 3>(arr2);
    view2.assign({
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
                  41,
                  42,
                  43,
                  44,
                  45,
                  });

    const auto printShape = [](auto& slice1, auto& slice2){
        Log::info() << slice1.rank() << std::endl;
        Log::info() << slice2.rank() << std::endl;
    };

    ArrayForEach<1>::apply(view1, view2, printShape);

}


} // namespace test
} // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
