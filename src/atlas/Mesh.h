/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Mesh_h
#define atlas_Mesh_h

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

#include "eckit/container/DenseMap.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/Metadata.h"
#include "atlas/Config.h"
#include "atlas/util/ObjectRegistry.h"

//------------------------------------------------------------------------------------------------------

// Forward declarations
namespace atlas { class Grid; }
namespace atlas { namespace mesh { class Nodes; } }
namespace atlas { namespace mesh { class HybridElements; } }
namespace atlas { class FunctionSpace; }
namespace atlas { class GridDistribution; }
namespace atlas { namespace mpl { class HaloExchange; } }
namespace atlas { namespace mpl { class GatherScatter; } }
namespace atlas { namespace mpl { class Checksum; } }

//----------------------------------------------------------------------------------------------------------------------


namespace atlas {

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
    FunctionSpace& create_function_space(const std::string& name,
                                         const std::string& shape_func,
                                         const std::vector<size_t>& shape);

    /// accessor by name
    FunctionSpace& function_space(const std::string& name) const;

    /// accessor by index
    FunctionSpace& function_space( size_t ) const;

    /// number of functional spaces
    size_t nb_function_spaces() const;
    
  protected:
    typedef eckit::DenseMap< std::string, eckit::SharedPtr<FunctionSpace> > StoreFS_t;

    StoreFS_t function_spaces_;

  };
}

class Mesh : public deprecated::FunctionSpaceContainer, 
             public util::Registered<Mesh> {

public: // types

    typedef eckit::SharedPtr<Mesh> Ptr;

public: // methods

    static Mesh* create( const eckit::Parametrisation& = Config() );
    static Mesh* create( const Grid&, const eckit::Parametrisation& = Config() );

    /// @brief Construct a empty Mesh
    explicit Mesh(const eckit::Parametrisation& = Config());

    /// @brief Construct mesh from grid.
    /// The mesh is global and only has a "nodes" FunctionSpace
    Mesh(const Grid&, const eckit::Parametrisation& = Config());

    /// Destructor
    /// @note No need to be virtual since this is not a base class.
    ~Mesh();

    Metadata& metadata() { return metadata_; }
    const Metadata& metadata() const { return metadata_; }

    void prettyPrint(std::ostream&) const;

    void print(std::ostream&) const;

    mesh::Nodes& createNodes(const Grid& g);

    const mesh::Nodes& nodes() const { return *nodes_; }
          mesh::Nodes& nodes()       { return *nodes_; }

    const mesh::HybridElements& cells() const { return *cells_; }
          mesh::HybridElements& cells()       { return *cells_; }

    const mesh::HybridElements& facets() const { return *facets_; }
          mesh::HybridElements& facets()       { return *facets_; }

    const mesh::HybridElements& ridges() const { return *ridges_; }
          mesh::HybridElements& ridges()       { return *ridges_; }

    const mesh::HybridElements& peaks() const { return *peaks_; }
          mesh::HybridElements& peaks()       { return *peaks_; }

    const mesh::HybridElements& edges() const { return *edges_; }
          mesh::HybridElements& edges()       { return *edges_; }

private:  // methods

    friend std::ostream& operator<<(std::ostream& s, const Mesh& p) {
        p.print(s);
        return s;
    }
    
    void createElements();

private: // members

    Metadata   metadata_;
    eckit::SharedPtr<mesh::Nodes> nodes_;           
                                                      // dimensionality : 2D    3D
    eckit::SharedPtr<mesh::HybridElements> cells_;    //                  2D    3D
    eckit::SharedPtr<mesh::HybridElements> facets_;   //                  1D    2D
    eckit::SharedPtr<mesh::HybridElements> ridges_;   //                  0D    1D
    eckit::SharedPtr<mesh::HybridElements> peaks_;    //                  NA    0D

    eckit::SharedPtr<mesh::HybridElements> edges_;  // alias to facets of 2D mesh, ridges of 3D mesh
    
    size_t dimensionality_;

public: // members to be removed
    void convert_old_to_new();
    bool has_grid() const { return grid_; }
    void set_grid( const Grid& p ) { grid_ = &p; }
    const Grid& grid() const {  ASSERT( grid_ ); return *grid_; }
private: // members to be removed
    const Grid* grid_;

};

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
#define mesh_Nodes mesh::Nodes
extern "C"
{
	Mesh* atlas__Mesh__new ();
	void atlas__Mesh__delete (Mesh* This);
  mesh_Nodes* atlas__Mesh__create_nodes (Mesh* This, int nb_nodes);
  void atlas__Mesh__create_function_space (Mesh* This, char* name,char* shape_func,int shape[], int shape_size, int fortran_ordering);
	FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name);
  mesh_Nodes* atlas__Mesh__nodes (Mesh* This);
}
#undef mesh_Nodes

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // atlas_Mesh_h
