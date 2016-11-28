/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestConnectivity
#include "ecbuild/boost_test_framework.h"
#include "atlas/mesh/Connectivity.h"

using namespace atlas::mesh;

namespace atlas {
namespace mesh {

#ifdef ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define INDEX_REF Index
#define FROM_FORTRAN -1
#define TO_FORTRAN   +1
#else
#define FORTRAN_BASE 0
#define INDEX_REF *
#define FROM_FORTRAN
#define TO_FORTRAN
#endif


BOOST_AUTO_TEST_CASE( test_irregular_connectivity )
{
    IrregularConnectivity conn("mesh");
    BOOST_CHECK_EQUAL(conn.rows(),0);
    BOOST_CHECK_EQUAL(conn.maxcols(),0);

    constexpr idx_t vals[4] = {2,3,5,6};
    conn.add(1, 4, vals, false);

    BOOST_CHECK_EQUAL(conn.rows(),1);
    BOOST_CHECK_EQUAL(conn.cols(0),4);
    BOOST_CHECK_EQUAL(conn.mincols(),4);
    BOOST_CHECK_EQUAL(conn.maxcols(),4);

    BOOST_CHECK_EQUAL(conn(0,0),3 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(0,1),4 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(0,2),6 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(0,3),7 + FROM_FORTRAN);

    BOOST_CHECK_EQUAL(conn.row(0)(0),3);
    BOOST_CHECK_EQUAL(conn.row(0)(1),4);
    BOOST_CHECK_EQUAL(conn.row(0)(2),6);
    BOOST_CHECK_EQUAL(conn.row(0)(3),7);


    constexpr idx_t vals2[6] = {1,3,4,3,7,8};
    conn.add(2, 3, vals2, true);

    BOOST_CHECK_EQUAL(conn.rows(),3);
    BOOST_CHECK_EQUAL(conn.cols(1),3);
    BOOST_CHECK_EQUAL(conn.cols(2),3);
    BOOST_CHECK_EQUAL(conn.mincols(),3);
    BOOST_CHECK_EQUAL(conn.maxcols(),4);

    BOOST_CHECK_EQUAL(conn(1,0),1 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(1,1),3 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(1,2),4 + FROM_FORTRAN);

    BOOST_CHECK_EQUAL(conn.row(2)(0),3);
    BOOST_CHECK_EQUAL(conn.row(2)(1),7);
    BOOST_CHECK_EQUAL(conn.row(2)(2),8);

    //set always add a TO_FORTRAN value
    conn.set(1,1,9);
    BOOST_CHECK_EQUAL(conn(1,0),1 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(1,1),9 + FROM_FORTRAN + TO_FORTRAN);
    BOOST_CHECK_EQUAL(conn(1,2),4 + FROM_FORTRAN);

    //set always add a TO_FORTRAN value
    constexpr idx_t vals3[3] = {6,7,5};
    conn.set(2, vals3);
    BOOST_CHECK_EQUAL(conn(2,0),6 + FROM_FORTRAN + TO_FORTRAN);
    BOOST_CHECK_EQUAL(conn(2,1),7 + FROM_FORTRAN + TO_FORTRAN);
    BOOST_CHECK_EQUAL(conn(2,2),5 + FROM_FORTRAN + TO_FORTRAN);


    constexpr idx_t vals4[8] = {2,11,51,12,4,13,55,78};
    conn.insert(1, 2, 4, vals4, false);
    BOOST_CHECK_EQUAL(conn.mincols(),3);
    BOOST_CHECK_EQUAL(conn.maxcols(),4);

    BOOST_CHECK_EQUAL(conn.rows(),5);
    BOOST_CHECK_EQUAL(conn.cols(1),4);
    BOOST_CHECK_EQUAL(conn.cols(2),4);

    BOOST_CHECK_EQUAL(conn(1,0),2 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(1,1),11 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(1,2),51 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(1,3),12 + FROM_FORTRAN );

    BOOST_CHECK_EQUAL(conn(2,0),4 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(2,1),13 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(2,2),55 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(2,3),78 + FROM_FORTRAN );

    BOOST_CHECK_EQUAL(conn(3,0),1 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(3,1),9 + FROM_FORTRAN + TO_FORTRAN);
    BOOST_CHECK_EQUAL(conn(3,2),4 + FROM_FORTRAN);

    constexpr idx_t vals5[2] = {3,6};
    conn.insert(3, 1, 2, vals5, true);

    BOOST_CHECK_EQUAL(conn.mincols(),2);
    BOOST_CHECK_EQUAL(conn.maxcols(),4);

    BOOST_CHECK_EQUAL(conn.rows(),6);
    BOOST_CHECK_EQUAL(conn.cols(3),2);
    BOOST_CHECK_EQUAL(conn.cols(4),3);

    BOOST_CHECK_EQUAL(conn(3,0),3 + FROM_FORTRAN + FORTRAN_BASE);
    BOOST_CHECK_EQUAL(conn(3,1),6 + FROM_FORTRAN + FORTRAN_BASE);

    BOOST_CHECK_EQUAL(conn(4,0),1 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(4,1),9 + FROM_FORTRAN + TO_FORTRAN);
    BOOST_CHECK_EQUAL(conn(4,2),4 + FROM_FORTRAN);

    //insert 3 rows with 1 column
    conn.insert(4, 3, 1);

    BOOST_CHECK_EQUAL(conn.rows(),9);
    BOOST_CHECK_EQUAL(conn.cols(4),1);
    BOOST_CHECK_EQUAL(conn.cols(5),1);
    BOOST_CHECK_EQUAL(conn.mincols(),1);
    BOOST_CHECK_EQUAL(conn.maxcols(),4);


    BOOST_CHECK_EQUAL(conn(7,0),1 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(7,1),9 + FROM_FORTRAN + TO_FORTRAN);
    BOOST_CHECK_EQUAL(conn(7,2),4 + FROM_FORTRAN);


    constexpr size_t cols[3] = {3,7,1};
    BOOST_CHECK_EQUAL(conn.cols(2),4);
    //insert in position 2, 3 rows with cols[3] number of columns
    conn.insert(2, 3, cols);

    BOOST_CHECK_EQUAL(conn.mincols(),1);
    BOOST_CHECK_EQUAL(conn.maxcols(),7);

    BOOST_CHECK_EQUAL(conn.rows(),12);
    BOOST_CHECK_EQUAL(conn.cols(2),3);
    BOOST_CHECK_EQUAL(conn.cols(3),7);
    BOOST_CHECK_EQUAL(conn.cols(4),1);
    BOOST_CHECK_EQUAL(conn.cols(5),4);

    BOOST_CHECK_EQUAL(conn(5,0),4 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(5,1),13 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(5,2),55 + FROM_FORTRAN );
    BOOST_CHECK_EQUAL(conn(5,3),78 + FROM_FORTRAN );

}

BOOST_AUTO_TEST_CASE( test_block_connectivity )
{
    idx_t vals[15] = {3,7,1,4,5,6,4,56,8,4,1,3,76,4,3};
    BlockConnectivity conn(3,5, vals);
    BOOST_CHECK_EQUAL(conn.rows(),3);
    BOOST_CHECK_EQUAL(conn.cols(),5);

    BOOST_CHECK_EQUAL(conn(0,2),1 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(1,1),4 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(2,2),76 + FROM_FORTRAN);

    idx_t vals2[10] = {2,3,9,34,356,86,3,24,84,45};

    conn.add(2,5, vals2);
    BOOST_CHECK_EQUAL(conn.rows(),5);
    BOOST_CHECK_EQUAL(conn.cols(),5);

    BOOST_CHECK_EQUAL(conn(0,2),1 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(1,1),4 + FROM_FORTRAN);
    BOOST_CHECK_EQUAL(conn(2,2),76 + FROM_FORTRAN);

    BOOST_CHECK_EQUAL(conn(3,2),9 + FROM_FORTRAN+FORTRAN_BASE);
    BOOST_CHECK_EQUAL(conn(3,4),356 + FROM_FORTRAN+FORTRAN_BASE);
    BOOST_CHECK_EQUAL(conn(4,1),3 + FROM_FORTRAN+FORTRAN_BASE);

}
BOOST_AUTO_TEST_CASE( test_block_connectivity_empty_add )
{
    BlockConnectivity conn;
    BOOST_CHECK_EQUAL(conn.rows(),0);
    BOOST_CHECK_EQUAL(conn.cols(),0);

    idx_t vals2[12] = {2,3,9,34,356,86,3,24,84,45,2,2};

    conn.add(2,5, vals2);
    BOOST_CHECK_EQUAL(conn.rows(),2);
    BOOST_CHECK_EQUAL(conn.cols(),5);

    BOOST_CHECK_EQUAL(conn(0,2),9 + FROM_FORTRAN+FORTRAN_BASE);
    BOOST_CHECK_EQUAL(conn(0,4),356 + FROM_FORTRAN+FORTRAN_BASE);
    BOOST_CHECK_EQUAL(conn(1,1),3 + FROM_FORTRAN+FORTRAN_BASE);
}


}
}
