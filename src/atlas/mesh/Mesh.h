/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>

#include "eckit/memory/SharedPtr.h"
#include "atlas/mesh/detail/MeshImpl.h"

//----------------------------------------------------------------------------------------------------------------------
// Forward declarations

namespace atlas {
    class Projection;
}

namespace atlas {
namespace util {
    class Metadata;
} }

namespace atlas {
namespace mesh {
    class Nodes;
    class HybridElements;
    typedef HybridElements Edges;
    typedef HybridElements Cells;
} }

namespace atlas {
namespace meshgenerator {
    class MeshGeneratorImpl;
} }

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

class Mesh {

public:

    using Implementation = mesh::detail::MeshImpl;
    using Nodes = mesh::Nodes;
    using Cells = mesh::Cells;
    using Edges = mesh::Edges;
    using HybridElements = mesh::HybridElements;
    using PartitionGraph = mesh::detail::PartitionGraph;
    using Polygon = mesh::PartitionPolygon;

public:

    Mesh();
    Mesh( const Mesh& );
    Mesh( const Implementation* );

    /// @brief Construct a mesh from a Stream (serialization)
    explicit Mesh(eckit::Stream&);

    /// @brief Serialization to Stream
    void encode(eckit::Stream& s) const { return impl_->encode(s); }

    void print(std::ostream& out) const { impl_->print(out); }

    const util::Metadata& metadata() const { return impl_->metadata(); }
          util::Metadata& metadata()       { return impl_->metadata(); }

    const Nodes& nodes() const { return impl_->nodes(); }
          Nodes& nodes()       { return impl_->nodes(); }

    const Cells& cells() const { return impl_->cells(); }
          Cells& cells()       { return impl_->cells();; }

    const Edges& edges() const { return impl_->edges(); }
          Edges& edges()       { return impl_->edges(); }

    const HybridElements& facets() const { return impl_->facets(); }
          HybridElements& facets()       { return impl_->facets(); }

    const HybridElements& ridges() const { return impl_->ridges(); }
          HybridElements& ridges()       { return impl_->ridges(); }

    const HybridElements& peaks() const { return impl_->peaks(); }
          HybridElements& peaks()       { return impl_->peaks(); }

    bool generated() const { return impl_->generated(); }

    /// @brief Return the memory footprint of the mesh
    size_t footprint() const { return impl_->footprint(); }

    size_t partition() const { return impl_->partition(); }

    size_t nb_partitions() const { return impl_->nb_partitions(); }

    void cloneToDevice() const { impl_->cloneToDevice(); }

    void cloneFromDevice() const { impl_->cloneFromDevice(); }

    void syncHostDevice() const { impl_->syncHostDevice(); }

    const Projection& projection() const { return impl_->projection(); }

    const PartitionGraph& partitionGraph() const { return impl_->partitionGraph(); }

    PartitionGraph::Neighbours nearestNeighbourPartitions() const { return impl_->nearestNeighbourPartitions(); }

    const Implementation* get() const { return impl_.get(); }
          Implementation* get()       { return impl_.get(); }

    const Polygon& polygon( size_t halo = 0) const { return impl_->polygon(halo); }

    const Grid& grid() const { return impl_->grid(); }

private:  // methods

    friend std::ostream& operator<<(std::ostream& s, const Mesh& p) {
        p.print(s);
        return s;
    }

    friend class meshgenerator::MeshGeneratorImpl;
    void setProjection(const Projection& p) { impl_->setProjection(p); }
    void setGrid(const Grid& p) { impl_->setGrid(p); }

private:

    eckit::SharedPtr<Implementation> impl_;

};

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas
