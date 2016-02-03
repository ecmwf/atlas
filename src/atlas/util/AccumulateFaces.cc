/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/FunctionSpace.h"
#include "atlas/util/AccumulateFaces.h"
#include "atlas/util/IndexView.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Nodes.h"

#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace util {

void accumulate_faces(
    FunctionSpace& func_space,
    std::vector< std::vector<int> >& node_to_face,
    std::vector<int>& face_nodes_data,
    std::vector< Face >& connectivity_edge_to_elem,
    size_t& nb_faces,
    size_t& nb_inner_faces )
{
  IndexView<int,2> elem_nodes( func_space.field( "nodes" ) );

    size_t nb_elems = func_space.shape(0);
    size_t nb_nodes_in_face = 2;

  std::vector< std::vector<int> > face_node_numbering;
  size_t nb_faces_in_elem;
  if (func_space.name() == "quads")
  {
    nb_faces_in_elem = 4;
    face_node_numbering.resize(nb_faces_in_elem, std::vector<int>(nb_nodes_in_face) );
    face_node_numbering[0][0] = 0;
    face_node_numbering[0][1] = 1;
    face_node_numbering[1][0] = 1;
    face_node_numbering[1][1] = 2;
    face_node_numbering[2][0] = 2;
    face_node_numbering[2][1] = 3;
    face_node_numbering[3][0] = 3;
    face_node_numbering[3][1] = 0;
  }
  else if (func_space.name() == "triags")
  {
    nb_faces_in_elem = 3;
    face_node_numbering.resize(nb_faces_in_elem, std::vector<int>(nb_nodes_in_face) );
    face_node_numbering[0][0] = 0;
    face_node_numbering[0][1] = 1;
    face_node_numbering[1][0] = 1;
    face_node_numbering[1][1] = 2;
    face_node_numbering[2][0] = 2;
    face_node_numbering[2][1] = 0;
  }
  else
  {
    throw eckit::BadParameter(func_space.name()+" is not \"quads\" or \"triags\"",Here());
  }

  for (size_t e = 0; e < nb_elems; ++e)
  {
    for (size_t f=0; f<nb_faces_in_elem; ++f)
    {
      bool found_face = false;

      std::vector<int> face_nodes(nb_nodes_in_face);
            for(size_t  jnode=0; jnode<nb_nodes_in_face; ++jnode)
        face_nodes[jnode] = elem_nodes(e,face_node_numbering[f][jnode]);

      int node = face_nodes[0];
            for(size_t jface=0; jface< node_to_face[node].size(); ++jface )
      {
        int face = node_to_face[node][jface];
        size_t nb_matched_nodes = 0;
        if (nb_nodes_in_face>1) // 2D or 3D
        {
          for(size_t jnode = 0; jnode < nb_nodes_in_face; ++jnode)
          {
                        size_t other_node = face_nodes[jnode];
                        for(size_t iface=0; iface<node_to_face[other_node].size(); ++iface )
            {
              if( node_to_face[face_nodes[jnode]][iface] == face )
              {
                ++nb_matched_nodes;
                break;
              }
            }
          }
          if (nb_matched_nodes == nb_nodes_in_face)
          {
            connectivity_edge_to_elem[face][1].f = func_space.index();
            connectivity_edge_to_elem[face][1].e = e;
            ++nb_inner_faces;
            found_face = true;
            break;
          }
        }
      }

      if (found_face == false)
      {
        connectivity_edge_to_elem.push_back( Face() );
        connectivity_edge_to_elem[nb_faces][0].f = func_space.index();
        connectivity_edge_to_elem[nb_faces][0].e = e;
        // if 2nd element stays negative, it is a bdry face
        connectivity_edge_to_elem[nb_faces][1].f = -1;
        connectivity_edge_to_elem[nb_faces][1].e = -1;
        for (size_t n = 0; n < nb_nodes_in_face; ++n)
        {
          node_to_face[face_nodes[n]].push_back(nb_faces);
          face_nodes_data.push_back(face_nodes[n]);
        }
        ++nb_faces;
      }
    }
  }
}

void accumulate_facets(
    const mesh::HybridElements &cells,
    const mesh::Nodes &nodes,
    std::vector< idx_t > &facet_nodes_data, // shape(nb_facets,nb_nodes_per_facet)
    std::vector< idx_t > &connectivity_facet_to_elem,
    size_t &nb_facets,
    size_t &nb_inner_facets,
    idx_t &missing_value )
{
  missing_value = -1;
  std::vector< std::vector<idx_t> > node_to_facet(nodes.size());
  nb_facets=0;
  nb_inner_facets=0;
  for( size_t t=0; t<cells.nb_types(); ++t )
  {
    const mesh::Elements& elements = cells.elements(t);
    const mesh::Elements::Connectivity& elem_nodes = elements.node_connectivity();

    size_t nb_elems = elements.size();
    size_t nb_nodes_in_facet = 2;

    std::vector< std::vector<int> > facet_node_numbering;
    size_t nb_facets_in_elem;
    if (elements.name() == "Quadrilateral")
    {
      nb_facets_in_elem = 4;
      facet_node_numbering.resize(nb_facets_in_elem, std::vector<int>(nb_nodes_in_facet) );
      facet_node_numbering[0][0] = 0;
      facet_node_numbering[0][1] = 1;
      facet_node_numbering[1][0] = 1;
      facet_node_numbering[1][1] = 2;
      facet_node_numbering[2][0] = 2;
      facet_node_numbering[2][1] = 3;
      facet_node_numbering[3][0] = 3;
      facet_node_numbering[3][1] = 0;
    }
    else if (elements.name() == "Triangle")
    {
      nb_facets_in_elem = 3;
      facet_node_numbering.resize(nb_facets_in_elem, std::vector<int>(nb_nodes_in_facet) );
      facet_node_numbering[0][0] = 0;
      facet_node_numbering[0][1] = 1;
      facet_node_numbering[1][0] = 1;
      facet_node_numbering[1][1] = 2;
      facet_node_numbering[2][0] = 2;
      facet_node_numbering[2][1] = 0;
    }
    else
    {
      throw eckit::BadParameter(elements.name()+" is not \"Quadrilateral\" or \"Triangle\"",Here());
    }

    std::vector<idx_t> facet_nodes(nb_nodes_in_facet);

    for (size_t e = 0; e < nb_elems; ++e)
    {
      for (size_t f=0; f<nb_facets_in_elem; ++f)
      {
        bool found_face = false;

        for(size_t  jnode=0; jnode<nb_nodes_in_facet; ++jnode){
          facet_nodes[jnode] = elem_nodes(e,facet_node_numbering[f][jnode]);
        }

        int node = facet_nodes[0];
        for(size_t jface=0; jface< node_to_facet[node].size(); ++jface )
        {
          int face = node_to_facet[node][jface];
          size_t nb_matched_nodes = 0;
          if (nb_nodes_in_facet>1) // 2D or 3D
          {
            for(size_t jnode = 0; jnode < nb_nodes_in_facet; ++jnode)
            {
              size_t other_node = facet_nodes[jnode];
              for(size_t iface=0; iface<node_to_facet[other_node].size(); ++iface )
              {
                if( node_to_facet[facet_nodes[jnode]][iface] == face )
                {
                  ++nb_matched_nodes;
                  break;
                }
              }
            }
            if (nb_matched_nodes == nb_nodes_in_facet)
            {
              connectivity_facet_to_elem[2*face+1] = e+elements.begin();
              ++nb_inner_facets;
              found_face = true;
              break;
            }
          }
        }

        if (found_face == false)
        {
          connectivity_facet_to_elem.push_back( elements.begin()+e     );
          // if 2nd element stays missing_value, it is a bdry face
          connectivity_facet_to_elem.push_back( missing_value );
          for (size_t n = 0; n < nb_nodes_in_facet; ++n)
          {
            node_to_facet[facet_nodes[n]].push_back(nb_facets);
            facet_nodes_data.push_back(facet_nodes[n]);
          }
          ++nb_facets;
        }
      }
    }
    }
}


} // namespace util
} // namespace atlas
