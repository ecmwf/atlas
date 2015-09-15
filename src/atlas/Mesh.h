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
namespace atlas { class Nodes; }
namespace atlas { class FunctionSpace; }
namespace atlas { class GridDistribution; }
namespace atlas { namespace mpl { class HaloExchange; } }
namespace atlas { namespace mpl { class GatherScatter; } }
namespace atlas { namespace mpl { class Checksum; } }

//----------------------------------------------------------------------------------------------------------------------


namespace atlas {

template <typename T>
class Store {
public:
  bool has(const std::string& name) const
  {
    return (store_.find(name) != store_.end());
  }
  void add( T* item )
  {
    ASSERT( !item->name().empty() );
    store_[item->name()] = eckit::SharedPtr<T>(item);
  }
  T& get(const std::string& name)
  {
    if( ! has(name) ) throw eckit::OutOfRange(name+" not found in store.",Here());
    return *store_.find(name)->second;
  }
  const T& get(const std::string& name) const
  {
    if( ! has(name) ) throw eckit::OutOfRange(name+" not found in store.",Here());
    return *store_.find(name)->second;
  }
  void remove(const std::string& name)
  {
    if( ! has(name) ) throw eckit::OutOfRange(name+" not found in store.",Here());
    store_.erase( store_.find(name) );
  }
private:
  std::map< std::string, eckit::SharedPtr<T> > store_;
};


class Mesh : public eckit::Owned, public util::Registered<Mesh> {

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

    /// checks if has a Grid
    bool has_grid() const { return grid_; }

    /// assign a Grid to this Mesh
    void set_grid( const Grid& p ) { grid_ = &p; }

    /// accessor of the Grid
    const Grid& grid() const {  ASSERT( grid_ ); return *grid_; }


    void prettyPrint(std::ostream&) const;

    void print(std::ostream&) const;


    Nodes& createNodes(const Grid& g);

    Nodes& createNodes( size_t );

    const Nodes& nodes() const { ASSERT(nodes_); return *nodes_; }
          Nodes& nodes()       { ASSERT(nodes_); return *nodes_; }

    const Store<const mpl::HaloExchange>& halo_exchange() const { return halo_exchange_; }
          Store<const mpl::HaloExchange>& halo_exchange()       { return halo_exchange_; }

    const Store<const mpl::GatherScatter>& gather_scatter() const { return gather_scatter_; }
          Store<const mpl::GatherScatter>& gather_scatter()       { return gather_scatter_; }

    const Store<const mpl::Checksum>& checksum() const { return checksum_; }
          Store<const mpl::Checksum>& checksum()       { return checksum_; }

private:  // methods

    friend std::ostream& operator<<(std::ostream& s, const Mesh& p) {
        p.print(s);
        return s;
    }

private: // members to be removed

    const Grid* grid_;

    typedef eckit::DenseMap< std::string, eckit::SharedPtr<FunctionSpace> > StoreFS_t;

    StoreFS_t function_spaces_;


private: // members

    Metadata   metadata_;
    eckit::SharedPtr<Nodes> nodes_;
    Store<const mpl::HaloExchange> halo_exchange_;
    Store<const mpl::GatherScatter> gather_scatter_;
    Store<const mpl::Checksum> checksum_;

};

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C"
{
	Mesh* atlas__Mesh__new ();
	void atlas__Mesh__delete (Mesh* This);
  Nodes* atlas__Mesh__create_nodes (Mesh* This, int nb_nodes);
  void atlas__Mesh__create_function_space (Mesh* This, char* name,char* shape_func,int shape[], int shape_size, int fortran_ordering);
	FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name);
  Nodes* atlas__Mesh__nodes (Mesh* This);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // atlas_Mesh_h
