/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_testing_AtlasFixture_h
#define atlas_testing_AtlasFixture_h

#include "eckit/testing/Setup.h"

#include "atlas/atlas.h"
#include "eckit/runtime/Main.h"
#include "eckit/mpi/Comm.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace test {


struct AtlasFixture : public eckit::testing::Setup {

    AtlasFixture()
    {
        eckit::Main::instance().taskID(
            eckit::mpi::comm("world").rank());
        if( eckit::Main::instance().taskID() != 0 )
            Log::reset();
        atlas_init();
    }

    ~AtlasFixture()
    {
        atlas_finalize();
    }
};

//----------------------------------------------------------------------------------------------------------------------

} // end namespace test
} // end namespace atlas

#endif // atlas_testing_AtlasFixture_h
