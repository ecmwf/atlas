/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Mesh_h
#define atlas_Mesh_h

#include <map>
#include <iosfwd>
#include <string>
#include <vector>

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/library/config.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Config.h"
#include "atlas/grid/Projection.h"
#include "atlas/parallel/mpi/mpi.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {
    class Grid;
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

namespace atlas {
namespace util {
namespace parallel {
namespace mpl {
    class HaloExchange;
    class GatherScatter;
    class Checksum;
} } } }

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace mesh {

class Mesh;
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

    void prettyPrint(std::ostream&) const;

    void print(std::ostream&) const;

    mesh::Nodes& createNodes(const grid::Grid& g);

    const mesh::Nodes& nodes() const { return *nodes_; }
          mesh::Nodes& nodes()       { return *nodes_; }

    const mesh::Cells& cells() const { return *cells_; }
          mesh::Cells& cells()       { return *cells_; }

    const mesh::Edges& edges() const { return *edges_; }
          mesh::Edges& edges()       { return *edges_; }

    const mesh::HybridElements& facets() const { return *facets_; }
          mesh::HybridElements& facets()       { return *facets_; }

    const mesh::HybridElements& ridges() const { return *ridges_; }
          mesh::HybridElements& ridges()       { return *ridges_; }

    const mesh::HybridElements& peaks() const { return *peaks_; }
          mesh::HybridElements& peaks()       { return *peaks_; }

    bool generated() const;

    /// @brief Return the memory footprint of the mesh
    size_t footprint() const;

    size_t nb_partitions() const;

    void cloneToDevice() const;

    void cloneFromDevice() const;

    void syncHostDevice() const;

    const grid::Projection& projection() const { return projection_; }

private:  // methods

    friend class Mesh;

    friend std::ostream& operator<<(std::ostream& s, const MeshImpl& p) {
        p.print(s);
        return s;
    }

    void createElements();

    friend class meshgenerator::MeshGeneratorImpl;
    void setProjection(const grid::Projection&);

private: // members

    util::Metadata   metadata_;

    eckit::SharedPtr<mesh::Nodes> nodes_;
                                                      // dimensionality : 2D | 3D
                                                      //                  --------
    eckit::SharedPtr<mesh::HybridElements> cells_;    //                  2D | 3D
    eckit::SharedPtr<mesh::HybridElements> facets_;   //                  1D | 2D
    eckit::SharedPtr<mesh::HybridElements> ridges_;   //                  0D | 1D
    eckit::SharedPtr<mesh::HybridElements> peaks_;    //                  NA | 0D

    eckit::SharedPtr<mesh::HybridElements> edges_;  // alias to facets of 2D mesh, ridges of 3D mesh

    size_t dimensionality_;

    grid::Projection projection_;
};

class Mesh {
public:

    using mesh_t = MeshImpl;

private:

    eckit::SharedPtr<MeshImpl> mesh_;

public:

    // operator MeshImpl&() { return *mesh_; }
    Mesh( const Mesh& other );
    Mesh( const MeshImpl* );
    Mesh();

    /// @brief Construct a mesh from a Stream (serialization)
    explicit Mesh(eckit::Stream&);

    /// @brief Serialization to Stream
    void encode(eckit::Stream& s) const { return mesh_->encode(s); }

    /// Destructor
    /// @note No need to be virtual since this is not a base class.
    ~Mesh() {}

    util::Metadata& metadata() { return mesh_->metadata(); }
    const util::Metadata& metadata() const { return mesh_->metadata(); }

    void prettyPrint(std::ostream& out) const { mesh_->prettyPrint(out); }

    void print(std::ostream& out) const { mesh_->print(out); }

    mesh::Nodes& createNodes(const grid::Grid& g) { return mesh_->createNodes(g); }

    const mesh::Nodes& nodes() const { return mesh_->nodes(); }
          mesh::Nodes& nodes()       { return mesh_->nodes(); }

    const mesh::Cells& cells() const { return mesh_->cells(); }
          mesh::Cells& cells()       { return mesh_->cells();; }

    const mesh::Edges& edges() const { return mesh_->edges(); }
          mesh::Edges& edges()       { return mesh_->edges(); }

    const mesh::HybridElements& facets() const { return mesh_->facets(); }
          mesh::HybridElements& facets()       { return mesh_->facets(); }

    const mesh::HybridElements& ridges() const { return mesh_->ridges(); }
          mesh::HybridElements& ridges()       { return mesh_->ridges(); }

    const mesh::HybridElements& peaks() const { return mesh_->peaks(); }
          mesh::HybridElements& peaks()       { return mesh_->peaks(); }

    bool generated() const { return mesh_->generated(); }

    /// @brief Return the memory footprint of the mesh
    size_t footprint() const { return mesh_->footprint(); }

    size_t nb_partitions() const { return mesh_->nb_partitions(); }

    void cloneToDevice() const { mesh_->cloneToDevice(); }

    void cloneFromDevice() const { mesh_->cloneFromDevice(); }

    void syncHostDevice() const { mesh_->syncHostDevice(); }

    const grid::Projection& projection() const { return mesh_->projection(); }

    mesh_t* get() { return mesh_.get(); }
    const mesh_t* get() const { return mesh_.get(); }

private:  // methods

    friend std::ostream& operator<<(std::ostream& s, const Mesh& p) {
        p.print(s);
        return s;
    }

    void createElements() { mesh_->createElements(); }

    friend class meshgenerator::MeshGeneratorImpl;
    void setProjection(const grid::Projection& p) { mesh_->setProjection(p); }

};

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C"
{
  Mesh::mesh_t* atlas__Mesh__new ();
  void atlas__Mesh__delete (Mesh::mesh_t* This);
  mesh::Nodes* atlas__Mesh__create_nodes (Mesh::mesh_t* This, int nb_nodes);
  mesh::Nodes* atlas__Mesh__nodes (Mesh::mesh_t* This);
  mesh::Edges* atlas__Mesh__edges (Mesh::mesh_t* This);
  mesh::Cells* atlas__Mesh__cells (Mesh::mesh_t* This);
  size_t atlas__Mesh__footprint (Mesh::mesh_t* This);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif // atlas_Mesh_h
