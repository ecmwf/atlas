/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "eckit/runtime/Tool.h"
#include "eckit/value/CompositeParams.h"

#include "atlas/array.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid.h"
#include "atlas/grid/Grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/ObjectHandle.h"

#include "tests/AtlasTestEnvironment.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::grid;

//-----------------------------------------------------------------------------

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

void take_array( const array::Array& arr ) {
    EXPECT( arr.size() == 20 );
}

class TakeArray {
public:
    TakeArray( const array::Array& arr ) { EXPECT( arr.size() == 20 ); }
};

//-----------------------------------------------------------------------------

CASE( "test_fieldcreator" ) {
    Field field( util::Config( "creator", "ArraySpec" )( "shape", array::make_shape( 10, 2 ) )(
        "datatype", array::DataType::real32().str() )( "name", "myfield" ) );

    EXPECT( field.datatype() == array::DataType::real32() );
    EXPECT( field.name() == "myfield" );

    Grid g( "O6" );

    Field arr( util::Config( "creator", "ArraySpec" )( "shape", array::make_shape( 10, 2 ) ) );
    EXPECT( arr.shape( 0 ) == 10 );
    EXPECT( arr.shape( 1 ) == 2 );
    EXPECT( arr.datatype() == array::DataType::real64() );

    util::Config ifs_parameters = util::Config( "creator", "IFS" )( "nlev", 137 )( "nproma", 10 )( "ngptot", g.size() );

    Log::info() << "Creating IFS field " << std::endl;
    Field ifs( util::Config( ifs_parameters )( "name", "myfield" )( "datatype", array::DataType::int32().str() )(
        "nvar", 8 ) );

    ATLAS_DEBUG_VAR( ifs );
    EXPECT( ifs.shape( 0 ) == 36 );
    EXPECT( ifs.shape( 1 ) == 8 );
    EXPECT( ifs.shape( 2 ) == 137 );
    EXPECT( ifs.shape( 3 ) == 10 );

    Log::flush();
}

CASE( "test_implicit_conversion" ) {
    Field field( "tmp", array::make_datatype<double>(), array::make_shape( 10, 2 ) );
    const array::Array& const_array = field;
    array::Array& array             = field;

    array::ArrayView<double, 2> arrv = array::make_view<double, 2>( array );
    arrv( 0, 0 )                     = 8.;

    array::ArrayView<const double, 2> carrv = array::make_view<double, 2>( const_array );
    EXPECT( carrv( 0, 0 ) == 8. );

    array::ArrayView<const double, 2> cfieldv = array::make_view<const double, 2>( field );
    EXPECT( cfieldv( 0, 0 ) == 8. );

    take_array( field );
    TakeArray ta( field );

    const Field& f = field;
    TakeArray cta( f );
}

CASE( "test_wrap_rawdata_through_array" ) {
    std::vector<double> rawdata( 20, 8. );
    util::ObjectHandle<array::Array> array( array::Array::wrap( rawdata.data(), array::make_shape( 10, 2 ) ) );
    Field field( "wrapped", array.get() );

    EXPECT( array->owners() == 2 );
    const array::ArrayView<double, 2> cfieldv = array::make_view<double, 2>( field );
    EXPECT( cfieldv( 9, 1 ) == 8. );
}

CASE( "test_wrap_rawdata_direct" ) {
    std::vector<double> rawdata( 10 * 2, 8. );
    Field field( "wrapped", rawdata.data(), array::make_shape( 10, 2 ) );

    EXPECT( field.array().owners() == 1 );
    const array::ArrayView<double, 2> cfieldv = array::make_view<double, 2>( field );
    EXPECT( cfieldv( 9, 1 ) == 8. );
}

CASE( "test_wrap_rawdata_through_field" ) {
    std::vector<double> rawdata( 10 * 2, 8. );
    Field field( "name", rawdata.data(), array::make_shape( 10, 2 ) );
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
