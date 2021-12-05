/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <memory>

#include "atlas/interpolation/method/Method.h"
#include "atlas/mesh/Halo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/KDTree.h"

namespace atlas {
namespace interpolation {
namespace method {

class KNearestNeighboursBase : public Method {
public:
    KNearestNeighboursBase(const Config& config): Method(config) {}
    virtual ~KNearestNeighboursBase() override {}

protected:
    void buildPointSearchTree(Mesh& meshSource) { buildPointSearchTree(meshSource, mesh::Halo(meshSource)); }
    void buildPointSearchTree(Mesh& meshSource, const mesh::Halo&);
    void buildPointSearchTree(const FunctionSpace&);
    bool extractTreeFromCache(const Cache&);

    util::IndexKDTree pTree_;
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
