/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/KNearestNeighboursBase.h"

#include "eckit/log/Timer.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/library/LibAtlas.h"


namespace atlas {
namespace interpolation {
namespace method {


void KNearestNeighboursBase::buildPointSearchTree(mesh::Mesh& meshSource) {
    using namespace atlas;
    eckit::TraceTimer<Atlas> tim("atlas::interpolation::method::KNearestNeighboursBase::setup()");


    // generate 3D point coordinates
    mesh::actions::BuildXYZField("xyz")(meshSource);
    array::ArrayView< double, 2 > coords(meshSource.nodes().field( "xyz" ));

    std::vector<PointIndex3::Value> pidx;
    pidx.reserve(meshSource.nodes().size());


    // build point-search tree
    pTree_.reset(new PointIndex3);

    for (size_t ip = 0; ip < meshSource.nodes().size(); ++ip) {
        PointIndex3::Point p(coords[ip].data());
        pidx.push_back(PointIndex3::Value(p, ip));
    }

    pTree_->build(pidx.begin(), pidx.end());
}


}  // method
}  // interpolation
}  // atlas
