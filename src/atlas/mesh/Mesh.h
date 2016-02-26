/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/atlas_config.h"
#include "atlas/internals/ObjectRegistry.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Config.h"

// -----------------------------------------------------------------------------
// Forward declarations
// -----------------------------------------------------------------------------
namespace atlas {
namespace grid {
    class Grid;
    class GridDistribution;
} }

namespace atlas {
namespace mesh {
    class Nodes;
    class HybridElements;
    typedef HybridElements Edges;
    typedef HybridElements Cells;
} }

namespace atlas {
namespace deprecated {
    class FunctionSpace;
} }

namespace atlas {
namespace util {
namespace parallel {
namespace mpl {
    class HaloExchange;
    class GatherScatter;
    class Checksum;
} } } }
// -----------------------------------------------------------------------------

namespace atlas {
namespace mesh {

#if !DEPRECATE_OLD_FUNCTIONSPACE
namespace deprecated {


  /**
   * @brief The FunctionSpaceContainer class
   * This class is a simple base class that will be removed soon, as
   * part of the new design.
   * Don't use any of these functions.
   */
  class FunctionSpaceContainer: public eckit::Owned {
  public:
    /// checks if function space exists
    bool has_function_space(const std::string& name) const;

    /// Takes ownership, and will be deleted automatically
    deprecated::FunctionSpace& create_function_space(const std::string& name,
                                         const std::string& shape_func,
                                         const std::vector<size_t>& shape);

    /// accessor by name
    deprecated::FunctionSpace& function_space(const std::string& name) const;

    /// accessor by index
    deprecated::FunctionSpace& function_space( size_t ) const;

    /// number of functional spaces
    size_t nb_function_spaces() const;

  protected:
    std::vector< eckit::SharedPtr<deprecated::FunctionSpace> >  function_spaces_;  ///< field handle storage
    std::map< std::string, size_t >                             index_;            ///< name-to-index map, to refer fields by name


  };
}

class Mesh : public deprecated::FunctionSpaceContainer,
             public internals::Registered<Mesh> {
#else
class Mesh : public eckit::Owned,
             public internals::Registered<Mesh> {
#endif

public: // types

    typedef eckit::SharedPtr<Mesh> Ptr;

public: // methods

    static Mesh* create( const eckit::Parametrisation& = util::Config() );
    static Mesh* create( const grid::Grid&, const eckit::Parametrisation& = util::Config() );

    /// @brief Construct a empty Mesh
    explicit Mesh(const eckit::Parametrisation& = util::Config());

    /// @brief Construct mesh from grid.
    /// The mesh is global and only has a "nodes" FunctionSpace
    Mesh(const grid::Grid&, const eckit::Parametrisation& = util::Config());

    /// Destructor
    /// @note No need to be virtual since this is not a base class.
    ~Mesh();

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


private:  // methods

    friend std::ostream& operator<<(std::ostream& s, const Mesh& p) {
        p.print(s);
        return s;
    }

    void createElements();

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

public: // members to be removed
#if ! DEPRECATE_OLD_FUNCTIONSPACE
    void convert_new_to_old();
#endif
    bool has_grid() const { return grid_; }
    void set_grid( const grid::Grid& p ) { grid_ = &p; }
    const grid::Grid& grid() const {  ASSERT( grid_ ); return *grid_; }
private: // members to be removed
    const grid::Grid* grid_;

};

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
#define mesh_Nodes mesh::Nodes
#define mesh_Edges mesh::Edges
#define mesh_Cells mesh::Cells
#define deprecated_FunctionSpace deprecated::FunctionSpace
extern "C"
{
  Mesh* atlas__Mesh__new ();
  void atlas__Mesh__delete (Mesh* This);
  mesh_Nodes* atlas__Mesh__create_nodes (Mesh* This, int nb_nodes);
  void atlas__Mesh__create_function_space (Mesh* This, char* name,char* shape_func,int shape[], int shape_size, int fortran_ordering);
  deprecated_FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name);
  mesh_Nodes* atlas__Mesh__nodes (Mesh* This);
  mesh_Edges* atlas__Mesh__edges (Mesh* This);
  mesh_Cells* atlas__Mesh__cells (Mesh* This);
}
#undef deprecated_FunctionSpace
#undef mesh_Nodes
#undef mesh_Edges
#undef mesh_Cells

//----------------------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif // atlas_Mesh_h
