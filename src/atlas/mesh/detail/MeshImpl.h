/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/util/Metadata.h"
#include "atlas/projection/Projection.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
  class Grid;
  class Mesh;
namespace mesh {
    class Nodes;
    class HybridElements;
    typedef HybridElements Edges;
    typedef HybridElements Cells;
} }

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace mesh {
namespace detail {

//----------------------------------------------------------------------------------------------------------------------

class PartitionGraph : public eckit::Owned {
public:
  using Neighbours = std::vector<size_t>;
public:
  PartitionGraph();
  PartitionGraph( size_t values[], size_t rows, size_t displs[], size_t counts[] );
  size_t footprint() const;
  size_t size() const;
  Neighbours nearestNeighbours(const size_t partition) const;
  size_t maximumNearestNeighbours() const;
  operator bool() const;
private:
  void print(std::ostream&) const;
  friend std::ostream& operator<<(std::ostream& s, const PartitionGraph& p) {
      p.print(s);
      return s;
  }
private:
  std::vector<size_t> counts_;
  std::vector<size_t> displs_;
  std::vector<size_t> values_;
  size_t maximum_nearest_neighbours_;
};

//----------------------------------------------------------------------------------------------------------------------

class MeshImpl : public eckit::Owned {

public: // methods

    /// @brief Construct a empty MeshImpl
    explicit MeshImpl();

    /// @brief Construct a mesh from a Stream (serialization)
    explicit MeshImpl(eckit::Stream&);

    /// @brief Serialization to Stream
    void encode(eckit::Stream& s) const;

    /// Destructor
    /// @note No need to be virtual since this is not a base class.
    ~MeshImpl();

    util::Metadata& metadata() { return metadata_; }
    const util::Metadata& metadata() const { return metadata_; }

    void print(std::ostream&) const;

    const Nodes& nodes() const { return *nodes_; }
          Nodes& nodes()       { return *nodes_; }

    const Cells& cells() const { return *cells_; }
          Cells& cells()       { return *cells_; }

    const Edges& edges() const { return *edges_; }
          Edges& edges()       { return *edges_; }

    const HybridElements& facets() const { return *facets_; }
          HybridElements& facets()       { return *facets_; }

    const HybridElements& ridges() const { return *ridges_; }
          HybridElements& ridges()       { return *ridges_; }

    const HybridElements& peaks() const { return *peaks_; }
          HybridElements& peaks()       { return *peaks_; }

    bool generated() const;

    /// @brief Return the memory footprint of the mesh
    size_t footprint() const;

    size_t partition() const;
    size_t nb_partitions() const;

    void cloneToDevice() const;

    void cloneFromDevice() const;

    void syncHostDevice() const;

    const Projection& projection() const { return projection_; }

    const PartitionGraph& partitionGraph() const;

    PartitionGraph::Neighbours nearestNeighbourPartitions() const;

private:  // methods

    friend class ::atlas::Mesh;

    friend std::ostream& operator<<(std::ostream& s, const MeshImpl& p) {
        p.print(s);
        return s;
    }

    void createElements();

    void setProjection(const Projection&);

    void build_partition_graph() const;

private: // members

    util::Metadata   metadata_;

    eckit::SharedPtr<Nodes> nodes_;
                                                      // dimensionality : 2D | 3D
                                                      //                  --------
    eckit::SharedPtr<HybridElements> cells_;    //                  2D | 3D
    eckit::SharedPtr<HybridElements> facets_;   //                  1D | 2D
    eckit::SharedPtr<HybridElements> ridges_;   //                  0D | 1D
    eckit::SharedPtr<HybridElements> peaks_;    //                  NA | 0D

    eckit::SharedPtr<HybridElements> edges_;  // alias to facets of 2D mesh, ridges of 3D mesh

    size_t dimensionality_;

    Projection projection_;

    mutable eckit::SharedPtr<PartitionGraph> partition_graph_;
};

//----------------------------------------------------------------------------------------------------------------------

} // namespace detail
} // namespace mesh
} // namespace atlas
