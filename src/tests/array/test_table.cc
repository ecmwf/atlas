/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestTable
#include "ecbuild/boost_test_framework.h"
#include "atlas/array/Table.h"
#include "atlas/runtime/Log.h"
#include "tests/AtlasFixture.h"


using namespace atlas::array;

namespace atlas {
namespace test {

#ifdef ATLAS_HAVE_FORTRAN
#define IN_FORTRAN -1
#else
#define IN_FORTRAN
#endif

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_table )
{
    Table table;
    BOOST_CHECK_EQUAL(table.rows(),0);
    BOOST_CHECK_EQUAL(table.maxcols(),0);

    constexpr idx_t vals[4] = {2,3,5,6};
    table.add(1, 4, vals, /* fortran-array = */ false);

    BOOST_CHECK_EQUAL(table.rows(),1);
    BOOST_CHECK_EQUAL(table.cols(0),4);
    BOOST_CHECK_EQUAL(table.mincols(),4);
    BOOST_CHECK_EQUAL(table.maxcols(),4);

    auto trv = make_table_view(table);
    BOOST_CHECK_EQUAL(trv(0,0),2);
    BOOST_CHECK_EQUAL(trv(0,1),3);
    BOOST_CHECK_EQUAL(trv(0,2),5);
    BOOST_CHECK_EQUAL(trv(0,3),6);

    BOOST_CHECK_EQUAL(trv.row(0)(0),2);
    BOOST_CHECK_EQUAL(trv.row(0)(1),3);
    BOOST_CHECK_EQUAL(trv.row(0)(2),5);
    BOOST_CHECK_EQUAL(trv.row(0)(3),6);

    constexpr idx_t vals2[6] = {1,3,4,3,7,8};
    table.add(2, 3, vals2, /* fortran-array = */ true);

    table.dump(Log::info()); Log::info() << std::endl;

    BOOST_CHECK_EQUAL(table.rows(),3);
    BOOST_CHECK_EQUAL(table.cols(1),3);
    BOOST_CHECK_EQUAL(table.cols(2),3);
    BOOST_CHECK_EQUAL(table.mincols(),3);
    BOOST_CHECK_EQUAL(table.maxcols(),4);

    trv = make_table_view(table);
    BOOST_CHECK_EQUAL(trv(1,0),1 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(1,1),3 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(1,2),4 IN_FORTRAN);

    BOOST_CHECK_EQUAL(trv.row(2)(0),3 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv.row(2)(1),7 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv.row(2)(2),8 IN_FORTRAN);

    auto twv = make_table_view<false>(table);
    twv(1,1) = 9 IN_FORTRAN;
    BOOST_CHECK_EQUAL(twv(1,0),1 IN_FORTRAN);
    BOOST_CHECK_EQUAL(twv(1,1),9 IN_FORTRAN);
    BOOST_CHECK_EQUAL(twv(1,2),4 IN_FORTRAN);

    constexpr idx_t vals3[3] = {6 IN_FORTRAN,7 IN_FORTRAN,5 IN_FORTRAN};
    twv.row(2) = vals3;
    BOOST_CHECK_EQUAL(twv(2,0),6 IN_FORTRAN);
    BOOST_CHECK_EQUAL(twv(2,1),7 IN_FORTRAN);
    BOOST_CHECK_EQUAL(twv(2,2),5 IN_FORTRAN);

    // TODO: following statement should not be allowed to compile!
    trv.row(2) = vals3;

    constexpr idx_t vals4[8] = {2,11,51,12,4,13,55,78};

    table.insert(1, 2, 4, vals4, /* fortran-array = */ false);

    BOOST_CHECK_EQUAL(table.mincols(),3);
    BOOST_CHECK_EQUAL(table.maxcols(),4);

    BOOST_CHECK_EQUAL(table.rows(),5);
    BOOST_CHECK_EQUAL(table.cols(0),4);
    BOOST_CHECK_EQUAL(table.cols(1),4);
    BOOST_CHECK_EQUAL(table.cols(2),4);
    BOOST_CHECK_EQUAL(table.cols(3),3);
    BOOST_CHECK_EQUAL(table.cols(4),3);


    trv = make_table_view(table);
    BOOST_CHECK_EQUAL(trv(1,0),2  );
    BOOST_CHECK_EQUAL(trv(1,1),11 );
    BOOST_CHECK_EQUAL(trv(1,2),51 );
    BOOST_CHECK_EQUAL(trv(1,3),12 );

    BOOST_CHECK_EQUAL(trv(2,0),4  );
    BOOST_CHECK_EQUAL(trv(2,1),13 );
    BOOST_CHECK_EQUAL(trv(2,2),55 );
    BOOST_CHECK_EQUAL(trv(2,3),78 );

    BOOST_CHECK_EQUAL(trv(3,0),1 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(3,1),9 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(3,2),4 IN_FORTRAN);


    constexpr idx_t vals5[2] = {3,6};
    table.insert(3, 1, 2, vals5, true);

    BOOST_CHECK_EQUAL(table.mincols(),2);
    BOOST_CHECK_EQUAL(table.maxcols(),4);

    BOOST_CHECK_EQUAL(table.rows(),6);
    BOOST_CHECK_EQUAL(table.cols(3),2);
    BOOST_CHECK_EQUAL(table.cols(4),3);

    trv = make_table_view(table);
    BOOST_CHECK_EQUAL(trv(3,0),3 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(3,1),6 IN_FORTRAN);

    BOOST_CHECK_EQUAL(trv(4,0),1 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(4,1),9 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(4,2),4 IN_FORTRAN);

    //insert 3 rows with 1 column
    table.insert(4, 3, 1);

    BOOST_CHECK_EQUAL(table.rows(),9);
    BOOST_CHECK_EQUAL(table.cols(4),1);
    BOOST_CHECK_EQUAL(table.cols(5),1);
    BOOST_CHECK_EQUAL(table.mincols(),1);
    BOOST_CHECK_EQUAL(table.maxcols(),4);

    trv = make_table_view(table);
    BOOST_CHECK_EQUAL(trv(7,0),1 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(7,1),9 IN_FORTRAN);
    BOOST_CHECK_EQUAL(trv(7,2),4 IN_FORTRAN);

    constexpr size_t cols[3] = {3,7,1};
    BOOST_CHECK_EQUAL(table.cols(2),4);
    //insert in position 2, 3 rows with cols[3] number of columns
    table.insert(2, 3, cols);

    BOOST_CHECK_EQUAL(table.mincols(),1);
    BOOST_CHECK_EQUAL(table.maxcols(),7);

    BOOST_CHECK_EQUAL(table.rows(),12);
    BOOST_CHECK_EQUAL(table.cols(2),3);
    BOOST_CHECK_EQUAL(table.cols(3),7);
    BOOST_CHECK_EQUAL(table.cols(4),1);
    BOOST_CHECK_EQUAL(table.cols(5),4);

    trv = make_table_view(table);
    BOOST_CHECK_EQUAL(trv(5,0),4  );
    BOOST_CHECK_EQUAL(trv(5,1),13 );
    BOOST_CHECK_EQUAL(trv(5,2),55 );
    BOOST_CHECK_EQUAL(trv(5,3),78 );

}

BOOST_AUTO_TEST_CASE( test_irregular_insert )
{
    Table table;
    BOOST_CHECK_EQUAL(table.rows(),0);
    BOOST_CHECK_EQUAL(table.maxcols(),0);

    constexpr idx_t vals[4] = {2,3,5,6};
    table.insert(0, 1, 4, vals, false);

    BOOST_CHECK_EQUAL(table.rows(),1);
    BOOST_CHECK_EQUAL(table.cols(0),4);
    BOOST_CHECK_EQUAL(table.mincols(),4);
    BOOST_CHECK_EQUAL(table.maxcols(),4);

    auto trv = make_table_view(table);
    BOOST_CHECK_EQUAL(trv(0,0),2 );
    BOOST_CHECK_EQUAL(trv(0,1),3 );
    BOOST_CHECK_EQUAL(trv(0,2),5 );
    BOOST_CHECK_EQUAL(trv(0,3),6 );

    BOOST_CHECK_EQUAL(trv.row(0)(0),2);
    BOOST_CHECK_EQUAL(trv.row(0)(1),3);
    BOOST_CHECK_EQUAL(trv.row(0)(2),5);
    BOOST_CHECK_EQUAL(trv.row(0)(3),6);
}

BOOST_AUTO_TEST_CASE( test_table_row )
{
    Table table;
    BOOST_CHECK_EQUAL(table.rows(),0);
    BOOST_CHECK_EQUAL(table.maxcols(),0);

    constexpr idx_t vals[4] = {2,3,5,6};
    table.insert(0, 1, 4, vals, /*fortran_array =*/false);

    auto twv = make_table_view<false>(table);
    auto trv = make_table_view<true >(table);

    auto row_write = twv.row(0);
    auto row_read  = trv.row(0);

    BOOST_CHECK_EQUAL( row_read(1),  3 );
    BOOST_CHECK_EQUAL( row_write(1), 3 );
    row_write(1) =  5;
    row_write(1) += 2;
    BOOST_CHECK_EQUAL( row_read(1),  7 );
    BOOST_CHECK_EQUAL( row_write(1), 7 );
}

} // namespace test
} // namespace atlas
