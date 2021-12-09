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

#include <string>

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/ObjectHandle.h"

namespace eckit {
class Hash;
class Parametrisation;
}  // namespace eckit

namespace atlas {
class Mesh;
class Grid;
class Projection;
}  // namespace atlas

namespace atlas {
namespace grid {
class Distribution;
class Partitioner;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace meshgenerator {
class MeshGeneratorImpl;
}

//----------------------------------------------------------------------------------------------------------------------

class MeshGenerator : DOXYGEN_HIDE(public util::ObjectHandle<meshgenerator::MeshGeneratorImpl>) {
public:
    using Parameters = atlas::util::Config;

public:
    using Handle::Handle;
    MeshGenerator(const std::string&, const eckit::Parametrisation& = util::NoConfig());
    MeshGenerator(const eckit::Parametrisation&);

    void hash(eckit::Hash&) const;

    Mesh generate(const Grid&, const grid::Distribution&) const;
    Mesh generate(const Grid&, const grid::Partitioner&) const;
    Mesh generate(const Grid&) const;

    Mesh operator()(const Grid&, const grid::Distribution&) const;
    Mesh operator()(const Grid&, const grid::Partitioner&) const;
    Mesh operator()(const Grid&) const;

    std::string type() const;
};

//----------------------------------------------------------------------------------------------------------------------

// Shorthands
class StructuredMeshGenerator : public MeshGenerator {
public:
    StructuredMeshGenerator(const eckit::Parametrisation& config = util::NoConfig()):
        MeshGenerator("structured", config) {}
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
