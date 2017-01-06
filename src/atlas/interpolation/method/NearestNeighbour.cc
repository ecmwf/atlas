/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/NearestNeighbour.h"

#include "eckit/log/Plural.h"
#include "eckit/log/Seconds.h"
#include "eckit/log/Timer.h"
#include "eckit/mpi/Comm.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/runtime/LibAtlas.h"
#include "atlas/runtime/Log.h"
#include "atlas/interpolation/method/PointIndex3.h"


namespace atlas {
namespace interpolation {
namespace method {


namespace {


MethodBuilder<NearestNeighbour> __builder("nearest-neighbour");


}  // (anonymous namespace)


void NearestNeighbour::execute(Method::Matrix& matrix, mesh::Mesh& meshSource, mesh::Mesh& meshTarget) const {
    using namespace atlas;
    eckit::TraceTimer<LibAtlas> timer("Nearest::execute");

    // generate 3D point coordinates
    mesh::actions::BuildXYZField("xyz")(meshSource);
    mesh::actions::BuildXYZField("xyz")(meshTarget);

    array::ArrayView< double, 2 > icoords(meshSource.nodes().field( "xyz" ));
    array::ArrayView< double, 2 > ocoords(meshTarget.nodes().field( "xyz" ));

    size_t inp_npts = meshSource.nodes().size();
    size_t out_npts = meshTarget.nodes().size();


    // build point-search tree
    eckit::ScopedPtr<PointIndex3> pTree(new PointIndex3);
    {
        std::vector<PointIndex3::Value> pidx;
        pidx.reserve(inp_npts);

        for (size_t ip = 0; ip < inp_npts; ++ip) {
            PointIndex3::Point p(icoords[ip].data());
            pidx.push_back(PointIndex3::Value(p, ip));
        }

        pTree->build(pidx.begin(), pidx.end());
    }


    // fill the sparse matrix
    std::vector< Triplet > weights_triplets;
    weights_triplets.reserve(out_npts);

    for (size_t ip = 0; ip < out_npts; ++ip) {

        if (ip && (ip % 1000 == 0)) {
            double rate = ip / timer.elapsed();
            Log::info() << eckit::BigNum(ip) << " (at " << rate << " points/s)..." << std::endl;
        }

        // find the closest input point to the output point
        PointIndex3::Point p(ocoords[ip].data());
        PointIndex3::NodeInfo nn = pTree->nearestNeighbour(p);
        size_t jp = nn.payload();

        // insert the weights into the interpolant matrix
        ASSERT(jp < inp_npts);
        weights_triplets.push_back(Triplet(ip, jp, 1));
    }

    // fill sparse matrix and return
    Matrix A(out_npts, inp_npts, weights_triplets);
    matrix.swap(A);
}


}  // method
}  // interpolation
}  // atlas
