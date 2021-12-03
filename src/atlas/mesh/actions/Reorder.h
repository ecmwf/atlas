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
#include <vector>

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/Factory.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

namespace atlas {
class Mesh;

namespace mesh {
namespace actions {

// ------------------------------------------------------------------

/// Base class for reordering the mesh nodes and elements
class ReorderImpl : public util::Object {
public:
    ReorderImpl() = default;

public:  // -- static functions --
    /// Reorder the nodes in the given mesh using a given order
    /// - All fields in mesh.nodes() are reordered.
    /// - mesh.cells().node_connectivity() gets updated
    /// - mesh.edges().node_connectivity() gets updated
    static void reorderNodes(Mesh& mesh, const std::vector<idx_t>& order);

    /// Reorder the cells by lowest node local index within each cell
    static void reorderCellsUsingNodes(Mesh& mesh);

    /// Reorder the edges by lowest node local index within each edge
    static void reorderEdgesUsingNodes(Mesh& mesh);

public:  // -- member functions --
    /// Reorder the nodes in the given mesh using the order computed with the computeNodesOrder function
    /// Then apply reorderCellsUsingNodes and reorderEdgesUsingNodes
    virtual void operator()(Mesh&);

    virtual std::vector<idx_t> computeNodesOrder(Mesh&) = 0;
};

//----------------------------------------------------------------------------------------------------------------------

class ReorderFactory : public util::Factory<ReorderFactory> {
public:
    static std::string className() { return "ReorderFactory"; }
    static const ReorderImpl* build(const eckit::Parametrisation&);
    using Factory::Factory;

private:
    virtual const ReorderImpl* make(const eckit::Parametrisation&) = 0;
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class ReorderBuilder : public ReorderFactory {
private:
    virtual const ReorderImpl* make(const eckit::Parametrisation& param) { return new T(param); }

public:
    using ReorderFactory::ReorderFactory;
};

//----------------------------------------------------------------------------------------------------------------------

class Reorder : DOXYGEN_HIDE(public util::ObjectHandle<ReorderImpl>) {
public:
    using Parameters = atlas::util::Config;

public:
    using Handle::Handle;
    Reorder(const eckit::Parametrisation& config): Handle(ReorderFactory::build(config)) {}

    /// Reorder the nodes in the given mesh using the order computed with the ReorderImpl::computeNodesOrder function
    /// Then apply ReorderImpl::reorderCellsUsingNodes and ReorderImpl::reorderEdgesUsingNodes
    void operator()(Mesh& mesh) { get()->operator()(mesh); }
};

//----------------------------------------------------------------------------------------------------------------------

/// Dummy implemenation of ReorderImpl which does nothing
class NoReorder : public ReorderImpl {
public:
    NoReorder(const eckit::Parametrisation&) {}
    void operator()(Mesh&) override {}
    std::vector<idx_t> computeNodesOrder(Mesh&) override { return {}; }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
