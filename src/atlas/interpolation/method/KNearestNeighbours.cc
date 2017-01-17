/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/method/KNearestNeighbours.h"

#include "eckit/log/Plural.h"
#include "eckit/log/Timer.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/runtime/LibAtlas.h"
#include "atlas/runtime/Log.h"


namespace atlas {
namespace interpolation {
namespace method {


namespace {


MethodBuilder<KNearestNeighbours> __builder("k-nearest-neighbours");


}  // (anonymous namespace)


KNearestNeighbours::KNearestNeighbours(const Method::Config& config) : KNearestNeighboursBase(config) {
    k_ = 1;
    config.get("k-nearest-neighbours", k_);
    ASSERT(k_);
}


void KNearestNeighbours::setup(mesh::Mesh& meshSource, mesh::Mesh& meshTarget) {
    using namespace atlas;


    // build point-search tree
    buildPointSearchTree(meshSource);
    ASSERT(pTree_);


    // generate 3D point coordinates
    mesh::actions::BuildXYZField("xyz")(meshTarget);
    array::ArrayView< double, 2 > coords(meshTarget.nodes().field( "xyz" ));

    size_t inp_npts = meshSource.nodes().size();
    size_t out_npts = meshTarget.nodes().size();


    // fill the sparse matrix
    std::vector< Triplet > weights_triplets;
    weights_triplets.reserve(out_npts * k_);
    {
        std::vector<double> weights;

        eckit::TraceTimer<LibAtlas> timer("atlas::interpolation::method::NearestNeighbour::setup()");
        for (size_t ip = 0; ip < out_npts; ++ip) {

            if (ip && (ip % 1000 == 0)) {
                double rate = ip / timer.elapsed();
                Log::debug() << eckit::BigNum(ip) << " (at " << rate << " points/s)..." << std::endl;
            }

            // find the closest input points to the output point
            PointIndex3::Point p(coords[ip].data());
            PointIndex3::NodeList nn = pTree_->kNearestNeighbours(p, k_);

            // calculate weights (individual and total, to normalise) using distance squared
            const size_t npts = nn.size();
            ASSERT(npts);
            weights.resize(npts, 0);

            double sum = 0;
            for (size_t j = 0; j < npts; ++j) {
                PointIndex3::Point np = nn[j].point();
                const double d2 = eckit::geometry::Point3::distance2(p, np);

                weights[j] = 1. / (1. + d2);
                sum += weights[j];
            }
            ASSERT(sum > 0);

            // insert weights into the matrix
            for (size_t j = 0; j < npts; ++j) {
                size_t jp = nn[j].payload();
                ASSERT(jp < inp_npts);
                weights_triplets.push_back(Triplet(ip, jp, weights[j]/sum));
            }
        }
    }


    // fill sparse matrix and return
    Matrix A(out_npts, inp_npts, weights_triplets);
    matrix_.swap(A);
}


}  // method
}  // interpolation
}  // atlas
