/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mpl/MPL.hpp"
#include "atlas/util/Array.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/mesh/Util.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"

using namespace atlas::meshgen;

namespace atlas {
namespace actions {


void distribute_mesh( Mesh& mesh )
{
  EqualAreaPartitioner partitioner( MPL::size() );
  int mypart = MPL::rank();

  FunctionSpace& nodes = mesh.function_space("nodes");
  int nb_nodes = nodes.extents()[0];
  ArrayView<double,2> latlon ( nodes.field("coordinates") );
  ArrayView<int,1>    part   ( nodes.field("partition")   );

  /*
  Create structure which we can partition with multiple keys (lat and lon)
  */
  std::vector<NodeInt> nodes_int(nb_nodes);

  int n=0;
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    nodes_int[n].x = microdeg(latlon(jnode,XX));
    nodes_int[n].y = microdeg(latlon(jnode,XX));
    nodes_int[n].n = n;
    ++n;
  }
  /*
  Assign partition to nodes
  */
  partitioner.partition(nb_nodes,nodes_int.data(),part.data());
  std::vector<NodeInt>().swap(nodes_int); // Deallocate completely


  int nb_keep_nodes = 0;
  std::vector<int> keep_nodes(nb_nodes,0);
  std::vector<int> node_loc(nb_nodes,-1);

  /*
  Find Elements that belong to mypart, and nodes that belong to that element
  */
  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata<int>("type") == Entity::ELEMS )
    {
      int nb_elems = elements.extents()[0];
      ArrayView<int,1> elem_gidx( elements.field("glb_idx") );
      ArrayView<int,1> elem_part( elements.field("partition") );
      IndexView<int,2> elem_nodes( elements.field("nodes") );
      int nb_nodes_per_elem = elem_nodes.extents()[1];
      int nb_keep_elems = 0;
      std::vector<int> keep_elems(nb_elems,0);
      std::vector<int> elem_loc(nb_elems,-1);
      for( int jelem = 0; jelem<nb_elems; ++jelem )
      {
        int maxpart = -1;
        for( int jnode=0; jnode<nb_nodes_per_elem; ++jnode )
        {
          maxpart = std::max(maxpart, part(elem_nodes(jelem,jnode)));
        }
        if( maxpart == mypart ) // keep element
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
      std::vector<int> new_elem_gidx(nb_keep_elems);
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
    }
  }
}

// ------------------------------------------------------------------

// C wrapper interfaces to C++ routines
void atlas__distribute_mesh (Mesh* mesh)
{
  distribute_mesh( *mesh );
}

// ------------------------------------------------------------------


} // namespace actions
} // namespace atlas
