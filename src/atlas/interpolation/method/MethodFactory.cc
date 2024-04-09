/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "MethodFactory.h"

// for static linking
#include "cubedsphere/CubedSphereBilinear.h"
#include "knn/GridBoxAverage.h"
#include "knn/GridBoxMaximum.h"
#include "knn/KNearestNeighbours.h"
#include "knn/NearestNeighbour.h"
#include "sphericalvector/SphericalVector.h"
#include "structured/Cubic2D.h"
#include "structured/Cubic3D.h"
#include "structured/Linear2D.h"
#include "structured/Linear3D.h"
#include "structured/QuasiCubic2D.h"
#include "structured/QuasiCubic3D.h"
#include "unstructured/FiniteElement.h"
#include "unstructured/UnstructuredBilinearLonLat.h"


namespace atlas {
namespace interpolation {

namespace {

void force_link() {
    static struct Link {
        Link() {
            MethodBuilder<method::UnstructuredBilinearLonLat>();
            MethodBuilder<method::FiniteElement>();
            MethodBuilder<method::KNearestNeighbours>();
            MethodBuilder<method::NearestNeighbour>();
            MethodBuilder<method::Linear2D>();
            MethodBuilder<method::Linear3D>();
            MethodBuilder<method::Cubic2D>();
            MethodBuilder<method::Cubic3D>();
            MethodBuilder<method::QuasiCubic2D>();
            MethodBuilder<method::QuasiCubic3D>();
            MethodBuilder<method::GridBoxAverage>();
            MethodBuilder<method::GridBoxMaximum>();
            MethodBuilder<method::CubedSphereBilinear>();
            MethodBuilder<method::SphericalVector>();
        }
    } link;
}

}  // namespace

Method* MethodFactory::build(const std::string& name, const Method::Config& config) {
    force_link();
    auto factory = get(name);
    return factory->make(config);
}

}  // namespace interpolation
}  // namespace atlas
