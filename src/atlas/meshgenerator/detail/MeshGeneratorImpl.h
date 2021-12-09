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

#include "atlas/util/Object.h"

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

//----------------------------------------------------------------------------------------------------------------------

class MeshGeneratorImpl : public util::Object {
public:
    MeshGeneratorImpl();

    virtual ~MeshGeneratorImpl();

    virtual void hash(eckit::Hash&) const = 0;

    virtual void generate(const Grid&, const grid::Partitioner&, Mesh&) const;
    virtual void generate(const Grid&, const grid::Distribution&, Mesh&) const = 0;
    virtual void generate(const Grid&, Mesh&) const                            = 0;

    Mesh generate(const Grid&, const grid::Partitioner&) const;
    Mesh generate(const Grid&, const grid::Distribution&) const;
    Mesh generate(const Grid&) const;

    Mesh operator()(const Grid&, const grid::Distribution&) const;
    Mesh operator()(const Grid&, const grid::Partitioner&) const;
    Mesh operator()(const Grid&) const;

    virtual std::string type() const = 0;

protected:
    void generateGlobalElementNumbering(Mesh& mesh) const;
    void setProjection(Mesh&, const Projection&) const;
    void setGrid(Mesh&, const Grid&, const grid::Distribution&) const;
    void setGrid(Mesh&, const Grid&, const std::string& distribution) const;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
