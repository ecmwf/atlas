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

#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"

namespace atlas {
class Mesh;
class Grid;
}  // namespace atlas

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class DelaunayMeshGenerator : public MeshGenerator::Implementation {
public:
    DelaunayMeshGenerator();
    DelaunayMeshGenerator(const eckit::Parametrisation& p);

    virtual ~DelaunayMeshGenerator() override;

    std::string type() const override { return "delaunay"; }

private:  // methods
    virtual void hash(eckit::Hash&) const override;

    using MeshGenerator::Implementation::generate;
    virtual void generate(const Grid&, const grid::Distribution&, Mesh&) const override;
    virtual void generate(const Grid&, Mesh&) const override;

    std::string mpi_comm_;
    int part_;
    bool remove_duplicate_points_;
    bool reshuffle_;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
