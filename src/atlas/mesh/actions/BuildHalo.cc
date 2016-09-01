/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


/// @warning Still doesn't know about periodic BC to enlarge Halo

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>

#include "eckit/memory/ScopedPtr.h"

#include "atlas/internals/atlas_config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/field/Field.h"
#include "atlas/internals/AccumulateFaces.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/internals/Parameters.h"
#include "atlas/internals/PeriodicTransform.h"
#include "atlas/internals/Unique.h"
#include "atlas/internals/LonLatMicroDeg.h"
#include "atlas/internals/Functions.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/Array.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"

using atlas::internals::accumulate_facets;
using atlas::internals::Topology;
using atlas::internals::PeriodicTransform;
using atlas::internals::UniqueLonLat;
using atlas::internals::LonLatMicroDeg;
using atlas::internals::microdeg;

namespace atlas {
namespace mesh {
namespace actions {


typedef gidx_t uid_t;

// ------------------------------------------------------------------
class BuildHaloHelper;

void increase_halo( Mesh& mesh );
void increase_halo_interior( BuildHaloHelper& );




class EastWest: public PeriodicTransform
{
public:
	EastWest()
	{
		x_translation_ = -360.;
	}
};


class WestEast: public PeriodicTransform
{
public:
	WestEast()
	{
		x_translation_ = 360.;
	}
};

typedef std::vector< std::vector<idx_t> > Node2Elem;

void build_lookup_node2elem( const Mesh& mesh, Node2Elem& node2elem )
{
  const mesh::Nodes& nodes  = mesh.nodes();

  node2elem.resize(nodes.size());
  for( size_t jnode=0; jnode<node2elem.size(); ++jnode )
    node2elem[jnode].clear();

  const mesh::HybridElements::Connectivity& elem_nodes = mesh.cells().node_connectivity();

  size_t nb_elems = mesh.cells().size();
  for (size_t elem=0; elem<nb_elems; ++elem)
  {
    for (size_t n=0; n<elem_nodes.cols(elem); ++n)
    {
      int node = elem_nodes(elem,n);
      node2elem[node].push_back( elem );
    }
  }
}


void accumulate_partition_bdry_nodes( Mesh& mesh, std::vector<int>& bdry_nodes )
{
  std::set<int> bdry_nodes_set;

  std::vector< idx_t > facet_nodes_data;
  std::vector< idx_t > connectivity_facet_to_elem;
  size_t nb_facets(0);
  size_t nb_inner_facets(0);
  idx_t missing_value;
  accumulate_facets(
      /*in*/  mesh.cells(),
      /*in*/  mesh.nodes(),
      /*out*/ facet_nodes_data, // shape(nb_facets,nb_nodes_per_facet)
      /*out*/ connectivity_facet_to_elem,
      /*out*/ nb_facets,
      /*out*/ nb_inner_facets,
      /*out*/ missing_value );

  array::ArrayView<idx_t,2> facet_nodes(facet_nodes_data.data(), array::make_shape(nb_facets,2));
  array::ArrayView<idx_t,2> facet_elem_connectivity(connectivity_facet_to_elem.data(), array::make_shape(nb_facets,2));

  for( size_t jface=0; jface<nb_facets; ++jface )
  {
    if( facet_elem_connectivity(jface,1) == missing_value )
    {
      for( size_t jnode=0; jnode<2; ++jnode) // 2 nodes per face
      {
        bdry_nodes_set.insert(facet_nodes(jface,jnode));
      }
    }
  }
  bdry_nodes = std::vector<int>( bdry_nodes_set.begin(), bdry_nodes_set.end());
}

template< typename Predicate >
std::vector<int> filter_nodes(std::vector<int> nodes, const Predicate& predicate )
{
  std::vector<int> filtered; filtered.reserve(nodes.size());
  for( size_t jnode=0; jnode<nodes.size(); ++jnode )
  {
    int inode = nodes[jnode];
    if( predicate(inode) )
      filtered.push_back(inode);
  }
  return filtered;
}

class Notification
{
public:
  void add_error(const std::string& note, const eckit::CodeLocation& loc )
  {
    notes.push_back( note + " @ " + std::string(loc) );
  }
  void add_error(const std::string& note )
  {
    notes.push_back( note );
  }

  bool error() const { return notes.size() > 0; }
  void reset() { notes.clear(); }

  std::string str() const
  {
    std::stringstream stream;
    for(size_t jnote = 0; jnote < notes.size(); ++jnote)
    {
      if ( jnote > 0 ) stream << "\n";
      stream << notes[jnote];
    }
    return stream.str();
  }

  operator std::string() const { return str(); }

private:
  friend std::ostream& operator<<(std::ostream& s, const Notification& notes) { s << notes.str();  return s; }

private:
  std::vector<std::string> notes;
};


typedef std::map<uid_t,int> Uid2Node;
void build_lookup_uid2node( Mesh& mesh, Uid2Node& uid2node )
{
  Notification notes;
  mesh::Nodes& nodes         = mesh.nodes();
  array::ArrayView<double,2> lonlat   ( nodes.lonlat() );
  array::ArrayView<gidx_t,1> glb_idx  ( nodes.global_index() );
  size_t nb_nodes = nodes.size();

  UniqueLonLat compute_uid(nodes);

  uid2node.clear();
  for( size_t jnode=0; jnode<nb_nodes; ++jnode )
  {
    uid_t uid = compute_uid(jnode);

    if( uid2node.count(uid) > 0 )
    {
      int other = uid2node[uid];
      std::stringstream msg;
      msg << "Node uid: " << uid << "   " << glb_idx(jnode)
          << " (" << lonlat(jnode,internals::LON) <<","<< lonlat(jnode,internals::LAT)<<")  has already been added as node "
          << glb_idx(other) << " (" << lonlat(other,internals::LON) <<","<< lonlat(other,internals::LAT)<<")";
      notes.add_error(msg.str());
    }
    uid2node[uid] = jnode;
  }
  if( notes.error() )
    throw eckit::SeriousBug(notes.str(),Here());
}

void accumulate_elements( const Mesh& mesh,
                          const array::ArrayView<uid_t,1>& node_uid,
                          const Uid2Node& uid2node,
                          const Node2Elem& node2elem,
                          std::vector<int>& found_elements,
                          std::set< uid_t >& new_nodes_uid )
{
  const mesh::HybridElements::Connectivity &elem_nodes = mesh.cells().node_connectivity();
  const array::ArrayView<int,1> elem_part ( mesh.cells().partition() );
  const array::ArrayView<gidx_t,1> ngidx ( mesh.nodes().global_index() );
  const array::ArrayView<gidx_t,1> egidx ( mesh.cells().global_index() );

  size_t nb_nodes = node_uid.size();

  std::set< int > found_elements_set;

  for( size_t jnode=0; jnode<nb_nodes; ++jnode )
  {
    uid_t uid = node_uid(jnode);

    int inode = -1;
    // search and get node index for uid
    Uid2Node::const_iterator found = uid2node.find(uid);
    if( found != uid2node.end() )
    {
      inode = found->second;
    }
    if( inode != -1 && size_t(inode) < node2elem.size() )
    {
      for(size_t jelem = 0; jelem < node2elem[inode].size(); ++jelem)
      {
        int e = node2elem[inode][jelem];
        if( size_t(elem_part(e)) == eckit::mpi::comm().rank() )
        {
          found_elements_set.insert( e );
        }
      }
    }
  }

  // found_bdry_elements_set now contains elements for the nodes
  found_elements = std::vector<int>( found_elements_set.begin(), found_elements_set.end());

  UniqueLonLat compute_uid(mesh.nodes());

  // Collect all nodes
  new_nodes_uid.clear();
  for(size_t jelem = 0; jelem < found_elements.size(); ++jelem)
  {
    size_t e = found_elements[jelem];
    size_t nb_elem_nodes = elem_nodes.cols(e);
    for( size_t n=0; n<nb_elem_nodes; ++n )
    {
      new_nodes_uid.insert( compute_uid(elem_nodes(e,n)) );
    }
  }

  // Remove nodes we already have
  for( size_t jnode=0; jnode<nb_nodes; ++jnode)
  {
    new_nodes_uid.erase( node_uid(jnode) ) ;
  }
}

class BuildHaloHelper
{
public:
  struct Buffers
  {
    std::vector< std::vector<int>    > node_part;

    std::vector< std::vector<int>    > node_ridx;

    std::vector< std::vector<int>    > node_flags;

    std::vector< std::vector<uid_t>  > node_glb_idx;

    std::vector< std::vector<double> > node_lonlat;

    std::vector< std::vector<uid_t>  > elem_glb_idx;

    std::vector< std::vector<uid_t>  > elem_nodes_id;

    std::vector< std::vector<int>    > elem_nodes_displs;

    std::vector< std::vector<int>    > elem_part;

    std::vector< std::vector<int>    > elem_type;

    Buffers(Mesh& mesh)
    {
      node_part.resize(eckit::mpi::comm().size());
      node_ridx.resize(eckit::mpi::comm().size());
      node_flags.resize(eckit::mpi::comm().size());
      node_glb_idx.resize(eckit::mpi::comm().size());
      node_lonlat.resize(eckit::mpi::comm().size());
      elem_glb_idx.resize(eckit::mpi::comm().size());
      elem_nodes_id.resize(eckit::mpi::comm().size());
      elem_nodes_displs.resize(eckit::mpi::comm().size());
      elem_part.resize(eckit::mpi::comm().size());
      elem_type.resize(eckit::mpi::comm().size());
    }

    void print( std::ostream& os )
    {
      os << "Nodes\n"
         << "-----\n";
      size_t n(0);
      for( size_t jpart=0; jpart<eckit::mpi::comm().size(); ++jpart )
      {
        for( size_t jnode=0; jnode<node_glb_idx[jpart].size(); ++jnode )
        {
          os << std::setw(4) << n++ << " : " << node_glb_idx[jpart][jnode] << "\n";
        }
      }
      os << std::flush;
      os << "Cells\n"
         << "-----\n";
      size_t e(0);
      for( size_t jpart=0; jpart<eckit::mpi::comm().size(); ++jpart )
      {
        for( size_t jelem=0; jelem<elem_glb_idx[jpart].size(); ++jelem )
        {
          os << std::setw(4) << e++ << " :  [ t" << elem_type[jpart][jelem] << " -- p" << elem_part[jpart][jelem] <<  "]  " << elem_glb_idx[jpart][jelem] << "\n";
        }
      }
      os << std::flush;
    }
  };

  static void all_to_all(Buffers& send, Buffers& recv)
  {
    eckit::mpi::comm().allToAll(send.node_glb_idx,  recv.node_glb_idx);
    eckit::mpi::comm().allToAll(send.node_part,     recv.node_part);
    eckit::mpi::comm().allToAll(send.node_ridx,     recv.node_ridx);
    eckit::mpi::comm().allToAll(send.node_flags,    recv.node_flags);
    eckit::mpi::comm().allToAll(send.node_lonlat,   recv.node_lonlat);
    eckit::mpi::comm().allToAll(send.elem_glb_idx,  recv.elem_glb_idx);
    eckit::mpi::comm().allToAll(send.elem_nodes_id, recv.elem_nodes_id);
    eckit::mpi::comm().allToAll(send.elem_part,     recv.elem_part);
    eckit::mpi::comm().allToAll(send.elem_type,     recv.elem_type);
    eckit::mpi::comm().allToAll(send.elem_nodes_displs, recv.elem_nodes_displs);
  }


public:
  Mesh& mesh;
  array::ArrayView<double,2> lonlat;
  array::ArrayView<gidx_t,1> glb_idx;
  array::ArrayView<int   ,1> part;
  array::IndexView<int   ,1> ridx;
  array::ArrayView<int   ,1> flags;
  array::ArrayView<int   ,1> ghost;
  mesh::HybridElements::Connectivity* elem_nodes;
  array::ArrayView<int,   1> elem_part;
  array::ArrayView<gidx_t,1> elem_glb_idx;

  std::vector<int> bdry_nodes;
  Node2Elem node_to_elem;
  Uid2Node uid2node;
  UniqueLonLat compute_uid;


public:
  BuildHaloHelper( Mesh& _mesh ): mesh(_mesh),
    compute_uid(mesh.nodes())
  {
    update();
  }

  void update()
  {
    compute_uid.update();
    mesh::Nodes& nodes         = mesh.nodes();
    lonlat   = array::ArrayView<double,2> ( nodes.lonlat() );
    glb_idx  = array::ArrayView<gidx_t,1> ( nodes.global_index() );
    part     = array::ArrayView<int   ,1> ( nodes.partition() );
    ridx     = array::IndexView<int   ,1> ( nodes.remote_index() );
    flags    = array::ArrayView<int   ,1> ( nodes.field("flags") );
    ghost    = array::ArrayView<int   ,1> ( nodes.ghost() );

    elem_nodes   = &mesh.cells().node_connectivity();
    elem_part    = array::ArrayView<int,1> ( mesh.cells().partition() );
    elem_glb_idx = array::ArrayView<gidx_t,1> ( mesh.cells().global_index() );
  }

  template< typename NodeContainer, typename ElementContainer >
  void fill_sendbuffer(Buffers& buf,const NodeContainer& nodes_uid, const ElementContainer& elems, const int p)
  {
    int nb_nodes = nodes_uid.size();
    buf.node_glb_idx[p].resize(nb_nodes);
    buf.node_part   [p].resize(nb_nodes);
    buf.node_ridx   [p].resize(nb_nodes);
    buf.node_flags  [p].resize(nb_nodes,Topology::NONE);
    buf.node_lonlat [p].resize(2*nb_nodes);

    int jnode=0;
    typename NodeContainer::iterator it;
    for( it=nodes_uid.begin(); it!=nodes_uid.end(); ++it, ++jnode )
    {
      uid_t uid = *it;

      Uid2Node::iterator found = uid2node.find( uid );
      if( found != uid2node.end() ) // Point exists inside domain
      {
        int node = found->second;
        buf.node_glb_idx[p][jnode]       = glb_idx(node);
        buf.node_part   [p][jnode]       = part   (node);
        buf.node_ridx   [p][jnode]       = ridx   (node);
        buf.node_lonlat [p][jnode*2+internals::LON] = lonlat (node,internals::LON);
        buf.node_lonlat [p][jnode*2+internals::LAT] = lonlat (node,internals::LAT);
        Topology::set(buf.node_flags[p][jnode],flags(node)|Topology::GHOST);
      }
      else
      {
        Log::warning() << "Node with uid " << uid << " needed by ["<<p<<"] was not found in ["<<eckit::mpi::comm().rank()<<"]." << std::endl;
        ASSERT(false);
      }
    }

    size_t nb_elems = elems.size();

    size_t nb_elem_nodes(0);
    for( size_t jelem=0; jelem<nb_elems; ++jelem )
    {
      size_t ielem = elems[jelem];
      nb_elem_nodes += elem_nodes->cols(ielem);
    }

    buf.elem_glb_idx [p].resize(nb_elems);
    buf.elem_part    [p].resize(nb_elems);
    buf.elem_type    [p].resize(nb_elems);
    buf.elem_nodes_id[p].resize(nb_elem_nodes);
    buf.elem_nodes_displs[p].resize(nb_elems);
    size_t jelemnode(0);
    for( size_t jelem=0; jelem<nb_elems; ++jelem )
    {
      buf.elem_nodes_displs[p][jelem] = jelemnode;
      size_t ielem = elems[jelem];
      buf.elem_glb_idx[p][jelem] = compute_uid( elem_nodes->row(ielem) );
      buf.elem_part   [p][jelem] = elem_part(ielem);
      buf.elem_type   [p][jelem] = mesh.cells().type_idx(ielem);
      for( size_t jnode=0; jnode<elem_nodes->cols(ielem); ++jnode )
        buf.elem_nodes_id[p][jelemnode++] = compute_uid( (*elem_nodes)(ielem,jnode) );
    }
  }

  template< typename NodeContainer, typename ElementContainer >
  void fill_sendbuffer(Buffers& buf,const NodeContainer& nodes_uid, const ElementContainer& elems, const PeriodicTransform& transform, int newflags, const int p)
  {
    int nb_nodes = nodes_uid.size();
    buf.node_glb_idx[p].resize(nb_nodes);
    buf.node_part   [p].resize(nb_nodes);
    buf.node_ridx   [p].resize(nb_nodes);
    buf.node_flags  [p].resize(nb_nodes,Topology::NONE);
    buf.node_lonlat [p].resize(2*nb_nodes);

    int jnode=0;
    typename NodeContainer::iterator it;
    for( it=nodes_uid.begin(); it!=nodes_uid.end(); ++it, ++jnode )
    {
      uid_t uid = *it;

      Uid2Node::iterator found = uid2node.find( uid );
      if( found != uid2node.end() ) // Point exists inside domain
      {
        int node = found->second;
        buf.node_part   [p][jnode]      = part   (node);
        buf.node_ridx   [p][jnode]      = ridx   (node);
        buf.node_lonlat [p][jnode*2+internals::LON] = lonlat (node,internals::LON);
        buf.node_lonlat [p][jnode*2+internals::LAT] = lonlat (node,internals::LAT);
        transform(&buf.node_lonlat[p][jnode*2],-1);
        // Global index of node is based on UID of destination
        buf.node_glb_idx[p][jnode]      = internals::unique_lonlat(&buf.node_lonlat [p][jnode*2]);
        Topology::set(buf.node_flags[p][jnode],newflags);
      }
      else
      {
        Log::warning() << "Node with uid " << uid << " needed by ["<<p<<"] was not found in ["<<eckit::mpi::comm().rank()<<"]." << std::endl;
        ASSERT(false);
      }
    }

    size_t nb_elems = elems.size();

    size_t nb_elem_nodes(0);
    for( size_t jelem=0; jelem<nb_elems; ++jelem )
    {
      size_t ielem = elems[jelem];
      nb_elem_nodes += elem_nodes->cols(ielem);
    }

    buf.elem_glb_idx [p].resize(nb_elems);
    buf.elem_part    [p].resize(nb_elems);
    buf.elem_type    [p].resize(nb_elems);
    buf.elem_nodes_id[p].resize(nb_elem_nodes);
    buf.elem_nodes_displs[p].resize(nb_elems);
    size_t jelemnode(0);
    for( size_t jelem=0; jelem<nb_elems; ++jelem )
    {
      buf.elem_nodes_displs[p][jelem] = jelemnode;
      size_t ielem = elems[jelem];
      buf.elem_part   [p][jelem] = elem_part(ielem);
      buf.elem_type   [p][jelem] = mesh.cells().type_idx(ielem);
      std::vector<double> crds(elem_nodes->cols(ielem)*2);
      for( int jnode=0; jnode<elem_nodes->cols(ielem); ++jnode)
      {
        double crd[] = { lonlat( (*elem_nodes)(ielem,jnode),internals::LON) , lonlat( (*elem_nodes)(ielem,jnode),internals::LAT) };
        transform(crd,-1);
        buf.elem_nodes_id[p][jelemnode++] = internals::unique_lonlat(crd);
        crds[jnode*2+internals::LON] = crd[internals::LON];
        crds[jnode*2+internals::LAT] = crd[internals::LAT];
      }
      // Global index of element is based on UID of destination
      buf.elem_glb_idx[p][jelem] = internals::unique_lonlat( crds.data(), elem_nodes->cols(ielem) );
    }

  }


  void add_nodes(Buffers& buf)
  {
    mesh::Nodes& nodes = mesh.nodes();
    int nb_nodes = nodes.size();

    // Nodes might be duplicated from different Tasks. We need to identify unique entries
    std::set<uid_t> node_uid;
    for( int jnode=0; jnode<nb_nodes; ++jnode )
    {
      node_uid.insert( compute_uid(jnode) );
    }
    std::vector< std::vector<int> > rfn_idx(eckit::mpi::comm().size());
    for(size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
    {
      rfn_idx[jpart].reserve(buf.node_glb_idx[jpart].size());
    }

    int nb_new_nodes=0;
    for(size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
    {
      for(size_t n = 0; n < buf.node_glb_idx[jpart].size(); ++n)
      {
        double crd[] = { buf.node_lonlat[jpart][n*2+internals::LON], buf.node_lonlat[jpart][n*2+internals::LAT] };
        bool inserted = node_uid.insert( internals::unique_lonlat(crd) ).second;
        if( inserted ) {
          rfn_idx[jpart].push_back(n);
        }
      }
      nb_new_nodes += rfn_idx[jpart].size();
    }

    // Resize nodes
    // ------------
    nodes.resize( nb_nodes+nb_new_nodes );
    flags   = array::ArrayView<int,   1>( nodes.field("flags") );
    glb_idx = array::ArrayView<gidx_t,1>( nodes.global_index() );
    part    = array::ArrayView<int,   1>( nodes.partition() );
    ridx    = array::IndexView<int,   1>( nodes.remote_index() );
    lonlat  = array::ArrayView<double,2>( nodes.lonlat() );
    ghost   = array::ArrayView<int,   1>( nodes.ghost() );

    compute_uid.update();

    // Add new nodes
    // -------------
    int new_node=0;
    for(size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
    {
      for(size_t n = 0; n < rfn_idx[jpart].size(); ++n)
      {
        int loc_idx = nb_nodes+new_node;
        Topology::reset(flags(loc_idx),buf.node_flags[jpart][rfn_idx[jpart][n]]);
        ghost  (loc_idx)    = Topology::check(flags(loc_idx),Topology::GHOST);
        glb_idx(loc_idx)    = buf.node_glb_idx [jpart][rfn_idx[jpart][n]];
        part   (loc_idx)    = buf.node_part    [jpart][rfn_idx[jpart][n]];
        ridx   (loc_idx)    = buf.node_ridx    [jpart][rfn_idx[jpart][n]];
        lonlat (loc_idx,internals::LON) = buf.node_lonlat  [jpart][rfn_idx[jpart][n]*2+internals::LON];
        lonlat (loc_idx,internals::LAT) = buf.node_lonlat  [jpart][rfn_idx[jpart][n]*2+internals::LAT];
        uid_t uid = compute_uid(loc_idx);

        // make sure new node was not already there
        Uid2Node::iterator found = uid2node.find(uid);
        if( found != uid2node.end() )
        {
          int other = found->second;
          std::stringstream msg;
          msg << "New node with uid " << uid << ":\n"  << glb_idx(loc_idx)
              << "("<<lonlat(loc_idx,internals::LON)<<","<<lonlat(loc_idx,internals::LAT)<<")\n";
          msg << "Existing already loc "<< other << "  :  " << glb_idx(other)
              << "("<<lonlat(other,internals::LON)<<","<<lonlat(other,internals::LAT)<<")\n";
          throw eckit::SeriousBug(msg.str(),Here());
        }
        uid2node[ uid ] = nb_nodes+new_node;
        ++new_node;
      }
    }
  }


  void add_elements(Buffers& buf)
  {
    // Elements might be duplicated from different Tasks. We need to identify unique entries
    std::set<uid_t> elem_uid;
    int nb_elems = mesh.cells().size();
    array::ArrayView<gidx_t,1> node_glb_idx ( mesh.nodes().global_index() );

    for( int jelem=0; jelem<nb_elems; ++jelem )
    {
      elem_uid.insert( compute_uid(elem_nodes->row(jelem)) );
    }

    std::vector< std::vector<int> > received_new_elems(eckit::mpi::comm().size());
    for(size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
    {
      received_new_elems[jpart].reserve(buf.elem_glb_idx[jpart].size());
    }

    size_t nb_new_elems(0);
    for(size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
    {
      for(size_t e = 0; e < buf.elem_glb_idx[jpart].size(); ++e)
      {
        bool inserted = elem_uid.insert( buf.elem_glb_idx[jpart][e] ).second;
        if( inserted )
        {
          received_new_elems[jpart].push_back(e);
        }
      }
      nb_new_elems += received_new_elems[jpart].size();
    }

    std::vector< std::vector< std::vector<int> > >
        elements_of_type( mesh.cells().nb_types(),
                          std::vector< std::vector<int> >( eckit::mpi::comm().size() ) );
    std::vector<size_t> nb_elements_of_type( mesh.cells().nb_types(), 0 );

    for(size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
    {
      for(size_t jelem = 0; jelem < received_new_elems[jpart].size(); ++jelem)
      {
        int ielem = received_new_elems[jpart][jelem];
        elements_of_type[buf.elem_type[jpart][ielem]][jpart] . push_back(ielem);
        ++nb_elements_of_type[buf.elem_type[jpart][ielem]];
      }
    }

    for( size_t t=0; t<mesh.cells().nb_types(); ++t )
    {
      const std::vector< std::vector<int> > &elems = elements_of_type[t];
      mesh::Elements& elements = mesh.cells().elements(t);

      // Add new elements
      mesh::Elements::Connectivity &node_connectivity = elements.node_connectivity();
      if( nb_elements_of_type[t] == 0 ) continue;
      size_t new_elems_pos = elements.add(nb_elements_of_type[t]);

      array::ArrayView<gidx_t,1> elem_type_glb_idx = elements.view<gidx_t,1>( mesh.cells().global_index() );
      array::ArrayView<int,   1> elem_type_part    = elements.view<int,1>( mesh.cells().partition() );
      array::ArrayView<int,   1> elem_type_halo    = elements.view<int,1>( mesh.cells().halo() );

      // Copy information in new elements
      size_t new_elem(0);
      for(size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
      {
        for(size_t e = 0; e < elems[jpart].size(); ++e)
        {
          size_t jelem = elems[jpart][e];
          elem_type_glb_idx(new_elems_pos+new_elem)   = buf.elem_glb_idx[jpart][jelem];
          elem_type_part   (new_elems_pos+new_elem)   = buf.elem_part[jpart][jelem];
          elem_type_halo   (new_elems_pos+new_elem)   = 1;
          for( size_t n=0; n<node_connectivity.cols(); ++n )
            node_connectivity.set(new_elems_pos+new_elem,n ,  uid2node[ buf.elem_nodes_id[jpart][ buf.elem_nodes_displs[jpart][jelem]+n] ] );
          ++new_elem;
        }
      }
    }
  }

  void add_buffers(Buffers& buf)
  {
    add_nodes(buf);
    add_elements(buf);
    update();
  }

};


void increase_halo_interior( BuildHaloHelper& helper )
{
  helper.update();
  if (helper.node_to_elem.size() == 0 )
    build_lookup_node2elem(helper.mesh,helper.node_to_elem);

  if( helper.uid2node.size() == 0 )
    build_lookup_uid2node(helper.mesh,helper.uid2node);


  // All buffers needed to move elements and nodes
  BuildHaloHelper::Buffers sendmesh(helper.mesh);
  BuildHaloHelper::Buffers recvmesh(helper.mesh);

  // 1) Find boundary nodes of this partition:

  accumulate_partition_bdry_nodes(helper.mesh,helper.bdry_nodes);
  const std::vector<int>& bdry_nodes = helper.bdry_nodes;

  // 2) Communicate uid of these boundary nodes to other partitions

  std::vector<uid_t> send_bdry_nodes_uid(bdry_nodes.size());
  for(size_t jnode = 0; jnode < bdry_nodes.size(); ++jnode)
    send_bdry_nodes_uid[jnode] = helper.compute_uid(bdry_nodes[jnode]);

  size_t size = eckit::mpi::comm().size();
  atlas::util::parallel::mpi::Buffer<uid_t,1> recv_bdry_nodes_uid_from_parts(size);

  eckit::mpi::comm().allGatherv(send_bdry_nodes_uid.begin(), send_bdry_nodes_uid.end(), recv_bdry_nodes_uid_from_parts);

  for (size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
  {
    // 3) Find elements and nodes completing these elements in
    //    other tasks that have my nodes through its UID

    array::ArrayView<uid_t,1> recv_bdry_nodes_uid = recv_bdry_nodes_uid_from_parts[jpart];

    std::vector<int>  found_bdry_elems;
    std::set< uid_t > found_bdry_nodes_uid;

    accumulate_elements(helper.mesh,recv_bdry_nodes_uid,
                        helper.uid2node,
                        helper.node_to_elem,
                        found_bdry_elems,
                        found_bdry_nodes_uid);

    // 4) Fill node and element buffers to send back
    helper.fill_sendbuffer(sendmesh, found_bdry_nodes_uid, found_bdry_elems, jpart);
  }

  // 5) Now communicate all buffers
  helper.all_to_all(sendmesh, recvmesh);

  // 6) Adapt mesh
  helper.add_buffers(recvmesh);
}

class PeriodicPoints {
public:
  PeriodicPoints(Mesh& mesh, int flag, size_t N)
  {
    flag_ = flag;
    N_ = N;
    flags_ = array::ArrayView<int,1> ( mesh.nodes().field("flags") );
  }

  bool operator()(int j) const
  {
    if( j>=N_ ) return false;
    if( Topology::check(flags_(j),flag_) ) return true;
    return false;
  }
private:
  int N_;
  int flag_;
  array::ArrayView<int,1> flags_;
};

void increase_halo_periodic( BuildHaloHelper& helper, const PeriodicPoints& periodic_points, const PeriodicTransform& transform, int newflags)
{
  helper.update();
  build_lookup_node2elem(helper.mesh,helper.node_to_elem);
  build_lookup_uid2node(helper.mesh,helper.uid2node);

  // All buffers needed to move elements and nodes
  BuildHaloHelper::Buffers sendmesh(helper.mesh);
  BuildHaloHelper::Buffers recvmesh(helper.mesh);

  // 1) Find boundary nodes of this partition:

  if( ! helper.bdry_nodes.size() )
    accumulate_partition_bdry_nodes(helper.mesh,helper.bdry_nodes);

  std::vector<int> bdry_nodes = filter_nodes(helper.bdry_nodes,periodic_points);

  // 2) Compute transformed uid of these boundary nodes and send to other partitions

  std::vector<uid_t> send_bdry_nodes_uid(bdry_nodes.size());
  for(size_t jnode = 0; jnode < bdry_nodes.size(); ++jnode)
  {
    double crd[] = { helper.lonlat(bdry_nodes[jnode],internals::LON), helper.lonlat(bdry_nodes[jnode],internals::LAT) };
    transform(crd,+1);
    send_bdry_nodes_uid[jnode] = internals::unique_lonlat(crd);
  }

  size_t size = eckit::mpi::comm().size();
  atlas::util::parallel::mpi::Buffer<uid_t,1> recv_bdry_nodes_uid_from_parts(size);

  eckit::mpi::comm().allGatherv(send_bdry_nodes_uid.begin(), send_bdry_nodes_uid.end(), recv_bdry_nodes_uid_from_parts);

  for (size_t jpart = 0; jpart < eckit::mpi::comm().size(); ++jpart)
  {
    // 3) Find elements and nodes completing these elements in
    //    other tasks that have my nodes through its UID

    array::ArrayView<uid_t,1> recv_bdry_nodes_uid = recv_bdry_nodes_uid_from_parts[jpart];

    std::vector<int>  found_bdry_elems;
    std::set< uid_t > found_bdry_nodes_uid;

    accumulate_elements(helper.mesh,recv_bdry_nodes_uid,
                        helper.uid2node,
                        helper.node_to_elem,
                        found_bdry_elems,
                        found_bdry_nodes_uid);

    // 4) Fill node and element buffers to send back
    helper.fill_sendbuffer(sendmesh, found_bdry_nodes_uid, found_bdry_elems, transform, newflags, jpart);
  }

  // 5) Now communicate all buffers
  helper.all_to_all(sendmesh, recvmesh);

  // 6) Adapt mesh

  helper.add_buffers(recvmesh);

}

void build_halo(Mesh& mesh, int nb_elems )
{
  int jhalo = 0;
  mesh.metadata().get("halo",jhalo);

  for( ; jhalo<nb_elems; ++jhalo )
  {
    size_t nb_nodes_before_halo_increase = mesh.nodes().size();

    BuildHaloHelper helper(mesh);

    increase_halo_interior( helper );

    PeriodicPoints westpts(mesh,Topology::PERIODIC|Topology::WEST,nb_nodes_before_halo_increase);

    increase_halo_periodic( helper, westpts, WestEast(), Topology::PERIODIC|Topology::WEST|Topology::GHOST );

    PeriodicPoints eastpts(mesh,Topology::PERIODIC|Topology::EAST,nb_nodes_before_halo_increase);

    increase_halo_periodic( helper, eastpts, EastWest(), Topology::PERIODIC|Topology::EAST|Topology::GHOST );

    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<jhalo+1<<"]";
    mesh.metadata().set(ss.str(),mesh.nodes().size());
  }

  mesh.metadata().set("halo",nb_elems);
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_halo ( Mesh* mesh, int nb_elems ) {
#undef ATLAS_ERROR_HANDLING
#define ATLAS_ERROR_HANDLING(x) x
  ATLAS_ERROR_HANDLING( build_halo(*mesh, nb_elems) );
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

