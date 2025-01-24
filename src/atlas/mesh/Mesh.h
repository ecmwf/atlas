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

#include <iosfwd>

#include "atlas/library/config.h"
#include "atlas/mesh/detail/MeshImpl.h"
#include "atlas/util/ObjectHandle.h"

//----------------------------------------------------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Projection;
class Grid;
}  // namespace atlas

namespace atlas {
namespace grid {
class Partitioner;
class Distribution;
}
}  // namespace atlas

namespace atlas {
namespace util {
class Metadata;
class PartitionPolygons;
}  // namespace util
}  // namespace atlas

namespace atlas {
namespace mesh {
class MeshBuilder;
class Nodes;
class HybridElements;
typedef HybridElements Edges;
typedef HybridElements Cells;
}  // namespace mesh
}  // namespace atlas

namespace atlas {
namespace meshgenerator {
class MeshGeneratorImpl;
}
}  // namespace atlas

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

class Mesh : DOXYGEN_HIDE(public util::ObjectHandle<mesh::detail::MeshImpl>) {
public:
    using Nodes          = mesh::Nodes;
    using Cells          = mesh::Cells;
    using Edges          = mesh::Edges;
    using HybridElements = mesh::HybridElements;
    using PartitionGraph = mesh::detail::PartitionGraph;
    using Polygon        = mesh::PartitionPolygon;
    using Polygons       = util::PartitionPolygons;

public:
    using Handle::Handle;
    Mesh();

    operator bool() const;

    /// @brief Generate a mesh from a Grid with recommended mesh generator and partitioner strategy
    Mesh(const Grid&, const eckit::Configuration& = util::NoConfig());

    Mesh(const Grid&, const grid::Partitioner&, const eckit::Configuration& = util::NoConfig());

    Mesh(const Grid&, const grid::Distribution&, const eckit::Configuration& = util::NoConfig());

    /// @brief Construct a mesh from a Stream (serialization)
    explicit Mesh(eckit::Stream&);

    /// @brief Serialization to Stream
    void encode(eckit::Stream& s) const { return get()->encode(s); }

    const util::Metadata& metadata() const { return get()->metadata(); }
    util::Metadata& metadata() { return get()->metadata(); }

    const Nodes& nodes() const { return get()->nodes(); }
    Nodes& nodes() { return get()->nodes(); }

    const Cells& cells() const { return get()->cells(); }
    Cells& cells() { return get()->cells(); }

    const Edges& edges() const { return get()->edges(); }
    Edges& edges() { return get()->edges(); }

    const HybridElements& facets() const { return get()->facets(); }
    HybridElements& facets() { return get()->facets(); }

    const HybridElements& ridges() const { return get()->ridges(); }
    HybridElements& ridges() { return get()->ridges(); }

    const HybridElements& peaks() const { return get()->peaks(); }
    HybridElements& peaks() { return get()->peaks(); }

    bool generated() const { return get()->generated(); }

    /// @brief Return the memory footprint of the mesh
    size_t footprint() const { return get()->footprint(); }

    idx_t part() const { return get()->part(); }

    idx_t nb_parts() const { return get()->nb_parts(); }

    [[deprecated("Use 'atlas::mesh::Mesh::nb_parts() instead")]] // added in v0.35.0
    idx_t nb_partitions() const { return nb_parts(); }

    std::string mpi_comm() const { return get()->mpi_comm(); }

    void updateDevice() const { get()->updateDevice(); }

    void updateHost() const { get()->updateHost(); }

    void syncHostDevice() const { get()->syncHostDevice(); }

    const Projection& projection() const { return get()->projection(); }

    const PartitionGraph& partitionGraph() const { return get()->partitionGraph(); }

    PartitionGraph::Neighbours nearestNeighbourPartitions() const { return get()->nearestNeighbourPartitions(); }

    const Polygon& polygon(idx_t halo = 0) const { return get()->polygon(halo); }
    const Polygons& polygons() const { return get()->polygons(); }

    const Grid grid() const { return get()->grid(); }

private:  // methods
    void print(std::ostream& out) const { get()->print(out); }

    friend std::ostream& operator<<(std::ostream& s, const Mesh& p) {
        p.print(s);
        return s;
    }

    friend class mesh::MeshBuilder;
    friend class meshgenerator::MeshGeneratorImpl;
    void setProjection(const Projection& p) { get()->setProjection(p); }
    void setGrid(const Grid& p) { get()->setGrid(p); }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
