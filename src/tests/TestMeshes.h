/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"

using namespace atlas;
using namespace atlas::grid;

namespace atlas {
namespace test {

Mesh generate_mesh(const StructuredGrid& grid) {
    auto config = util::Config("partitioner", "equal_regions");
    StructuredMeshGenerator generate(config);
    return generate(grid);
}

Mesh generate_mesh(std::initializer_list<long> nx) {
    return generate_mesh(ReducedGaussianGrid(nx));
}

}  // end namespace test
}  // end namespace atlas
