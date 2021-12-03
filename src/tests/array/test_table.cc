/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/Table.h"
#include "atlas/library/defines.h"
#include "atlas/runtime/Log.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

#if ATLAS_HAVE_FORTRAN
#define IN_FORTRAN -1
#else
#define IN_FORTRAN
#endif

CASE("test_table") {
    Table table;
    EXPECT(table.rows() == 0);
    EXPECT(table.maxcols() == 0);

    constexpr idx_t vals[4] = {2, 3, 5, 6};
    table.add(1, 4, vals, /* fortran-array = */ false);

    EXPECT(table.rows() == 1);
    EXPECT(table.cols(0) == 4);
    EXPECT(table.mincols() == 4);
    EXPECT(table.maxcols() == 4);

    auto trv = make_table_view(table);
    EXPECT(trv(0, 0) == 2);
    EXPECT(trv(0, 1) == 3);
    EXPECT(trv(0, 2) == 5);
    EXPECT(trv(0, 3) == 6);

    EXPECT(trv.row(0)(0) == 2);
    EXPECT(trv.row(0)(1) == 3);
    EXPECT(trv.row(0)(2) == 5);
    EXPECT(trv.row(0)(3) == 6);

    constexpr idx_t vals2[6] = {1, 3, 4, 3, 7, 8};
    table.add(2, 3, vals2, /* fortran-array = */ true);

    table.dump(Log::info());
    Log::info() << std::endl;

    EXPECT(table.rows() == 3);
    EXPECT(table.cols(1) == 3);
    EXPECT(table.cols(2) == 3);
    EXPECT(table.mincols() == 3);
    EXPECT(table.maxcols() == 4);

    trv = make_table_view(table);
    EXPECT(trv(1, 0) == 1 IN_FORTRAN);
    EXPECT(trv(1, 1) == 3 IN_FORTRAN);
    EXPECT(trv(1, 2) == 4 IN_FORTRAN);

    EXPECT(trv.row(2)(0) == 3 IN_FORTRAN);
    EXPECT(trv.row(2)(1) == 7 IN_FORTRAN);
    EXPECT(trv.row(2)(2) == 8 IN_FORTRAN);

    auto twv  = make_table_view<false>(table);
    twv(1, 1) = 9 IN_FORTRAN;
    EXPECT(twv(1, 0) == 1 IN_FORTRAN);
    EXPECT(twv(1, 1) == 9 IN_FORTRAN);
    EXPECT(twv(1, 2) == 4 IN_FORTRAN);

    constexpr idx_t vals3[3] = {6 IN_FORTRAN, 7 IN_FORTRAN, 5 IN_FORTRAN};
    twv.row(2)               = vals3;
    EXPECT(twv(2, 0) == 6 IN_FORTRAN);
    EXPECT(twv(2, 1) == 7 IN_FORTRAN);
    EXPECT(twv(2, 2) == 5 IN_FORTRAN);

    // TODO: following statement should not be allowed to compile!
    trv.row(2) = vals3;

    constexpr idx_t vals4[8] = {2, 11, 51, 12, 4, 13, 55, 78};

    table.insert(1, 2, 4, vals4, /* fortran-array = */ false);

    EXPECT(table.mincols() == 3);
    EXPECT(table.maxcols() == 4);

    EXPECT(table.rows() == 5);
    EXPECT(table.cols(0) == 4);
    EXPECT(table.cols(1) == 4);
    EXPECT(table.cols(2) == 4);
    EXPECT(table.cols(3) == 3);
    EXPECT(table.cols(4) == 3);

    trv = make_table_view(table);
    EXPECT(trv(1, 0) == 2);
    EXPECT(trv(1, 1) == 11);
    EXPECT(trv(1, 2) == 51);
    EXPECT(trv(1, 3) == 12);

    EXPECT(trv(2, 0) == 4);
    EXPECT(trv(2, 1) == 13);
    EXPECT(trv(2, 2) == 55);
    EXPECT(trv(2, 3) == 78);

    EXPECT(trv(3, 0) == 1 IN_FORTRAN);
    EXPECT(trv(3, 1) == 9 IN_FORTRAN);
    EXPECT(trv(3, 2) == 4 IN_FORTRAN);

    constexpr idx_t vals5[2] = {3, 6};
    table.insert(3, 1, 2, vals5, true);

    EXPECT(table.mincols() == 2);
    EXPECT(table.maxcols() == 4);

    EXPECT(table.rows() == 6);
    EXPECT(table.cols(3), 2);
    EXPECT(table.cols(4), 3);

    trv = make_table_view(table);
    EXPECT(trv(3, 0) == 3 IN_FORTRAN);
    EXPECT(trv(3, 1) == 6 IN_FORTRAN);

    EXPECT(trv(4, 0) == 1 IN_FORTRAN);
    EXPECT(trv(4, 1) == 9 IN_FORTRAN);
    EXPECT(trv(4, 2) == 4 IN_FORTRAN);

    // insert 3 rows with 1 column
    table.insert(4, 3, 1);

    EXPECT(table.rows() == 9);
    EXPECT(table.cols(4) == 1);
    EXPECT(table.cols(5) == 1);
    EXPECT(table.mincols() == 1);
    EXPECT(table.maxcols() == 4);

    trv = make_table_view(table);
    EXPECT(trv(7, 0) == 1 IN_FORTRAN);
    EXPECT(trv(7, 1) == 9 IN_FORTRAN);
    EXPECT(trv(7, 2) == 4 IN_FORTRAN);

    constexpr size_t cols[3] = {3, 7, 1};
    EXPECT(table.cols(2), 4);
    // insert in position 2, 3 rows with cols[3] number of columns
    table.insert(2, 3, cols);

    EXPECT(table.mincols() == 1);
    EXPECT(table.maxcols() == 7);

    EXPECT(table.rows() == 12);
    EXPECT(table.cols(2) == 3);
    EXPECT(table.cols(3) == 7);
    EXPECT(table.cols(4) == 1);
    EXPECT(table.cols(5) == 4);

    trv = make_table_view(table);
    EXPECT(trv(5, 0) == 4);
    EXPECT(trv(5, 1) == 13);
    EXPECT(trv(5, 2) == 55);
    EXPECT(trv(5, 3) == 78);
}

CASE("test_irregular_insert") {
    Table table;
    EXPECT(table.rows() == 0);
    EXPECT(table.maxcols() == 0);

    constexpr idx_t vals[4] = {2, 3, 5, 6};
    table.insert(0, 1, 4, vals, false);

    EXPECT(table.rows() == 1);
    EXPECT(table.cols(0) == 4);
    EXPECT(table.mincols() == 4);
    EXPECT(table.maxcols() == 4);

    auto trv = make_table_view(table);
    EXPECT(trv(0, 0) == 2);
    EXPECT(trv(0, 1) == 3);
    EXPECT(trv(0, 2) == 5);
    EXPECT(trv(0, 3) == 6);

    EXPECT(trv.row(0)(0) == 2);
    EXPECT(trv.row(0)(1) == 3);
    EXPECT(trv.row(0)(2) == 5);
    EXPECT(trv.row(0)(3) == 6);
}

CASE("test_table_row") {
    Table table;
    EXPECT(table.rows() == 0);
    EXPECT(table.maxcols() == 0);

    constexpr idx_t vals[4] = {2, 3, 5, 6};
    table.insert(0, 1, 4, vals, /*fortran_array =*/false);

    auto twv = make_table_view<false>(table);
    auto trv = make_table_view<true>(table);

    auto row_write = twv.row(0);
    auto row_read  = trv.row(0);

    EXPECT(row_read(1) == 3);
    EXPECT(row_write(1) == 3);
    row_write(1) = 5;
    row_write(1) += 2;
    EXPECT(row_read(1) == 7);
    EXPECT(row_write(1) == 7);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
