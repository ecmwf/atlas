/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <limits>
#include "atlas/mpi/mpi.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"
#include "atlas/Util.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/meshgen/EqualRegionsPartitioner.h"

using namespace atlas::meshgen;

namespace atlas {
namespace actions {


void distribute_mesh( Mesh& mesh )
{
  int mypart = eckit::mpi::rank();

  FunctionSpace& nodes = mesh.function_space("nodes");
  int nb_nodes = nodes.shape(0);
  ArrayView<double,2> lonlat    ( nodes.field("lonlat") );
  ArrayView<int,   1> node_part ( nodes.field("partition")   );
  ArrayView<gidx_t,1> node_gidx ( nodes.field("glb_idx")   );

  const Grid& g = mesh.grid();
  EqualRegionsPartitioner partitioner( g );
  partitioner.partition(node_part.data());

  int nb_keep_nodes = 0;
  std::vector<int> keep_nodes(nb_nodes,0);
  std::vector<int> node_loc(nb_nodes,-1);

  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    if( node_part(jnode) == mypart )
    {
      keep_nodes[jnode] = 1;
      node_loc[jnode] = nb_keep_nodes;
      ++nb_keep_nodes;
    }
  }


  /*
  Find Elements that belong to mypart, and nodes that belong to that element
  */
  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata().get<long>("type") == Entity::ELEMS )
    {
      int nb_elems = elements.shape(0);
      ArrayView<gidx_t,1> elem_gidx( elements.field("glb_idx") );
      ArrayView<int,1> elem_part( elements.field("partition") );
      IndexView<int,2> elem_nodes( elements.field("nodes") );
      int nb_nodes_per_elem = elem_nodes.shape(1);
      int nb_keep_elems = 0;
      std::vector<int> keep_elems(nb_elems,0);
      std::vector<int> elem_loc(nb_elems,-1);
      for( int jelem = 0; jelem<nb_elems; ++jelem )
      {
        int maxpart = -1;
        int minpart = std::numeric_limits<int>::max();
        double y(0.);
        for( int jnode=0; jnode<nb_nodes_per_elem; ++jnode )
        {
          int n = elem_nodes(jelem,jnode);
          maxpart = std::max(maxpart, node_part(n));
          minpart = std::min(minpart, node_part(n));
          y += lonlat(n,LAT);
        }
        y /= static_cast<double>(nb_nodes_per_elem);
        if( (y>=0 && maxpart == mypart) || (y<0 && minpart == mypart) ) // keep element
        {
          keep_elems[jelem] = 1;
          elem_loc[jelem] = nb_keep_elems;
          ++nb_keep_elems;
          for( int jnode=0; jnode<nb_nodes_per_elem; ++jnode )
          {
            int p = elem_nodes(jelem,jnode);
            if( keep_nodes[p] == 0 ) // keep node
            {
              keep_nodes[p] = 1;
              node_loc[p] = nb_keep_nodes;
              ++nb_keep_nodes;
            }
          }
        }
      }
      std::vector<gidx_t> new_elem_gidx(nb_keep_elems);
      Array<int> new_elem_nodes_arr(nb_keep_elems,nb_nodes_per_elem);
      ArrayView<int,2> new_elem_nodes( new_elem_nodes_arr );
      for( int jelem = 0; jelem<nb_elems; ++jelem )
      {
        if( keep_elems[jelem] == 1 )
        {
          int e = elem_loc[jelem];
          for( int jnode=0; jnode<nb_nodes_per_elem; ++jnode )
          {
            new_elem_nodes(e,jnode) = node_loc[elem_nodes(jelem,jnode)];
          }
          new_elem_gidx[e] = elem_gidx(jelem);
        }
      }
      nb_elems = nb_keep_elems;
      std::vector<int> shape = elements.shape();
      shape[0] = nb_elems;
      elements.resize(shape);
      elem_gidx  = ArrayView<gidx_t,1> ( elements.field("glb_idx") );
      elem_part  = ArrayView<int,1> ( elements.field("partition") );
      elem_nodes = IndexView<int,2> ( elements.field("nodes") );
      for( int jelem=0; jelem<nb_elems; ++jelem )
      {
        elem_gidx(jelem) = new_elem_gidx[jelem];
        elem_part(jelem) = mypart;
        for( int jnode=0; jnode<nb_nodes_per_elem; ++jnode )
        {
          elem_nodes(jelem,jnode) = new_elem_nodes(jelem,jnode);
        }
      }
      elements.metadata().set("nb_owned",nb_elems);
    }
  }
  std::vector<int> new_node_gidx(nb_keep_nodes);
  std::vector<int> new_node_part(nb_keep_nodes);
  Array<double> new_lonlat_arr(nb_keep_nodes,2);
  ArrayView<double,2> new_lonlat(new_lonlat_arr);
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    if( keep_nodes[jnode] == 1 )
    {
      int inode = node_loc[jnode];
      new_node_gidx[inode] = node_gidx(jnode);
      new_node_part[inode] = node_part(jnode);
      new_lonlat(inode,LON) = lonlat(jnode,LON);
      new_lonlat(inode,LAT) = lonlat(jnode,LAT);
    }
  }
  nb_nodes = nb_keep_nodes;
  std::vector<int> shape = nodes.shape();
  shape[0] = nb_nodes;
  nodes.resize(shape);
  node_gidx  = ArrayView<gidx_t,   1> ( nodes.field("glb_idx") );
  node_part  = ArrayView<int,   1> ( nodes.field("partition") );
  lonlat     = ArrayView<double,2> ( nodes.field("lonlat") );
  int nb_owned = 0;
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    node_gidx(jnode) = new_node_gidx[jnode];
    node_part(jnode) = new_node_part[jnode];
    lonlat(jnode,LON) = new_lonlat(jnode,LON);
    lonlat(jnode,LAT) = new_lonlat(jnode,LAT);
    if( node_part(jnode) == mypart ) ++nb_owned;
  }
  nodes.metadata().set("nb_owned",nb_owned);
}

// ------------------------------------------------------------------

// C wrapper interfaces to C++ routines
void atlas__distribute_mesh (Mesh* mesh)
{
  ATLAS_ERROR_HANDLING( distribute_mesh( *mesh ) );
}

// ------------------------------------------------------------------


} // namespace actions
} // namespace atlas
