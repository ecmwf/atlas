/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Trace.h"
#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

CASE( "test elapsed" ) {
    auto trace = Trace( Here() );

    EXPECT( trace.running() );

    EXPECT( trace.elapsed() == 0. );

    trace.pause();

    EXPECT( trace.running() );

    double elapsed = trace.elapsed();
    EXPECT( elapsed != 0. );
    EXPECT( trace.elapsed() == elapsed );

    trace.resume();

    EXPECT( trace.running() );

    trace.stop();

    EXPECT( trace.elapsed() != elapsed );
}

CASE( "test trace OpenMP" ) {
    atlas_omp_parallel_for( int i = 0; i < 10; ++i ) {
        auto trace = Trace( Here(), "loop" );
        if ( ATLAS_HAVE_OMP ) {
            trace.stop();
            if ( atlas_omp_get_thread_num() > 0 ) {
                EXPECT( trace.elapsed() == 0. );
            }
            else {
                EXPECT( trace.elapsed() != 0. );
            }
        }
    }
}

CASE( "test barrier" ) {
    EXPECT( runtime::trace::Barriers::state() == Library::instance().traceBarriers() );
    {
        runtime::trace::Barriers set_barriers( true );
        EXPECT( runtime::trace::Barriers::state() == true );
        {
            runtime::trace::Barriers set_barriers( false );
            EXPECT( runtime::trace::Barriers::state() == false );
        }
        EXPECT( runtime::trace::Barriers::state() == true );
    }
    EXPECT( runtime::trace::Barriers::state() == Library::instance().traceBarriers() );
}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
