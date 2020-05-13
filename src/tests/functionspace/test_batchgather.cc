/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/log/Bytes.h"
#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_functionspace_StructuredColumns_batchgather" ) {
    size_t root          = 0;
    std::string gridname = eckit::Resource<std::string> ("--grid", "O8");
    int nfields = eckit::Resource<int> ("--fields", 3);
    int irank = mpi::comm ().rank ();
    int nproc = mpi::comm ().size ();

    std::cout << " nfields = " << nfields << std::endl;
    Grid grid (gridname);


    Distribution dist (grid, grid::Partitioner ("equal_regions"));
    functionspace::StructuredColumns fs (grid, dist, util::Config ("halo", 1) 
                                       | util::Config ("periodic_points", true));

    std::vector<int> partition (grid.size ());

    for (int i = 0; i < grid.size (); i++)
      partition[i] = dist.partition (i);

    FieldSet sloc;
    FieldSet sglo;

    for (int i = 0; i < nfields; i++)
      {
        std::string name = std::string ("#") + std::to_string (i);
        Field floc = fs.createField<double> (option::name (name) | option::name ("owner", (i % nproc)));
        auto v = array::make_view<double,1> (floc);
        for (int j = 0; j < fs.sizeOwned (); j++)
          v[j] = irank + static_cast<double> (j) / 1000.0;
        sloc.add (floc);
        Field fglo = Field (name, atlas::array::DataType::kind<double> (), atlas::array::make_shape (grid.size ()));
        sglo.add (fglo);
      }

    fs.gather (sloc, sglo);

    for (int i = 0; i < nfields; i++)
      {
        const auto & fglo = sglo[i];
        const auto v = array::make_view<double,1> (fglo);
        int owner = -1;
        EXPECT (f.metadata ().get ("owner", owner));
        if (owner == irank)
          {
            for (int j = 0; j < grid.size (); j++)
              
          }
      }

}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
