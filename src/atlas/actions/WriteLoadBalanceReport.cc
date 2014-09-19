/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iomanip>
#include <fstream>
#include "eckit/filesystem/LocalPathName.h"
#include "atlas/mpl/MPL.h"
#include "atlas/mesh/Util.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/actions/WriteLoadBalanceReport.h"

namespace atlas {
namespace actions {

void write_load_balance_report( const Mesh& mesh, const std::string& filename )
{
  int npart = MPL::size();
  int root = 0;

  std::vector<int> nb_total_nodes(npart,0);
  std::vector<int> nb_owned_nodes(npart,0);
  std::vector<int> nb_ghost_nodes(npart,0);

  std::vector<int> nb_total_edges(npart,0);
  std::vector<int> nb_owned_edges(npart,0);
  std::vector<int> nb_ghost_edges(npart,0);

  if( mesh.has_function_space("nodes") )
  {
    FunctionSpace& nodes = mesh.function_space("nodes");
    IsGhost is_ghost(nodes);
    int nb_nodes = nodes.shape(0);
    int nowned(0);
    int nghost(0);
    for( int n=0; n<nb_nodes; ++n )
    {
      if( is_ghost(n) )
        ++nghost;
      else
        ++nowned;
    }
    MPL_CHECK_RESULT( MPI_Gather( &nb_nodes, 1, MPI_INT,
                                  nb_total_nodes.data(), 1, MPI_INT,
                                  root, MPI_COMM_WORLD ) );
    MPL_CHECK_RESULT( MPI_Gather( &nowned, 1, MPI_INT,
                                  nb_owned_nodes.data(), 1, MPI_INT,
                                  root, MPI_COMM_WORLD ) );
    MPL_CHECK_RESULT( MPI_Gather( &nghost, 1, MPI_INT,
                                  nb_ghost_nodes.data(), 1, MPI_INT,
                                  root, MPI_COMM_WORLD ) );  }

  if( mesh.has_function_space("edges") )
  {
    FunctionSpace& nodes = mesh.function_space("nodes");
    IsGhost is_ghost(nodes);
    FunctionSpace& edges = mesh.function_space("edges");
    IndexView<int,2> edge_nodes ( edges.field("nodes") );
    int nb_edges = edges.shape(0);
    int nowned(0);
    int nghost(0);
    for( int j=0; j<nb_edges; ++j )
    {
      if( is_ghost(edge_nodes(j,0)) )
        ++nghost;
      else
        ++nowned;
    }
    MPL_CHECK_RESULT( MPI_Gather( &nb_edges, 1, MPI_INT,
                                  nb_total_edges.data(), 1, MPI_INT,
                                  root, MPI_COMM_WORLD ) );
    MPL_CHECK_RESULT( MPI_Gather( &nowned, 1, MPI_INT,
                                  nb_owned_edges.data(), 1, MPI_INT,
                                  root, MPI_COMM_WORLD ) );
    MPL_CHECK_RESULT( MPI_Gather( &nghost, 1, MPI_INT,
                                  nb_ghost_edges.data(), 1, MPI_INT,
                                  root, MPI_COMM_WORLD ) );
  }

  if( MPL::rank() == 0 )
  {
    eckit::LocalPathName path(filename);
    std::ofstream ofs;
    int idt = 8;
    ofs.open( path.c_str(), std::ofstream::out );
    ofs << "# STATISTICS\n";
    ofs << std::setw(1)  << "#" << std::setw(5) << "";
    ofs << std::setw(idt) << "nodes";
    ofs << std::setw(idt) << "onodes";
    ofs << std::setw(idt) << "gnodes";
    ofs << std::setw(idt) << "edges";
    ofs << std::setw(idt) << "oedges";
    ofs << std::setw(idt) << "gedges";
    ofs << "\n";
    ofs << std::setw(6)  << "# tot";
    ofs << std::setw(idt) << std::accumulate(nb_total_nodes.data(),nb_total_nodes.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_owned_nodes.data(),nb_owned_nodes.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_ghost_nodes.data(),nb_ghost_nodes.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_total_edges.data(),nb_total_edges.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_owned_edges.data(),nb_owned_edges.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_ghost_edges.data(),nb_ghost_edges.data()+npart,0);
    ofs << "\n";
    ofs << std::setw(6)  << "# max";
    ofs << std::setw(idt) << *std::max_element(nb_total_nodes.data(),nb_total_nodes.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_owned_nodes.data(),nb_owned_nodes.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_ghost_nodes.data(),nb_ghost_nodes.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_total_edges.data(),nb_total_edges.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_owned_edges.data(),nb_owned_edges.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_ghost_edges.data(),nb_ghost_edges.data()+npart);
    ofs << "\n";
    ofs << std::setw(6)  << "# min";
    ofs << std::setw(idt) << *std::min_element(nb_total_nodes.data(),nb_total_nodes.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_owned_nodes.data(),nb_owned_nodes.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_ghost_nodes.data(),nb_ghost_nodes.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_total_edges.data(),nb_total_edges.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_owned_edges.data(),nb_owned_edges.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_ghost_edges.data(),nb_ghost_edges.data()+npart);
    ofs << "\n";
    ofs << std::setw(6)  << "# avg";
    ofs << std::setw(idt) << std::accumulate(nb_total_nodes.data(),nb_total_nodes.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_owned_nodes.data(),nb_owned_nodes.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_ghost_nodes.data(),nb_ghost_nodes.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_total_edges.data(),nb_total_edges.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_owned_edges.data(),nb_owned_edges.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_ghost_edges.data(),nb_ghost_edges.data()+npart,0)/npart;
    ofs << "\n";
    ofs << "#----------------------------------------------------\n";
    ofs << "# PER TASK\n";
    ofs << std::setw(6)  << "# part";
    ofs << std::setw(idt) << "nodes";
    ofs << std::setw(idt) << "onodes";
    ofs << std::setw(idt) << "gnodes";
    ofs << std::setw(idt) << "edges";
    ofs << std::setw(idt) << "oedges";
    ofs << std::setw(idt) << "gedges";
    ofs << "\n";
    for( int jpart=0; jpart<npart; ++jpart )
    {
      ofs << std::setw(6)  << jpart;
      ofs << std::setw(idt) << nb_total_nodes[jpart];
      ofs << std::setw(idt) << nb_owned_nodes[jpart];
      ofs << std::setw(idt) << nb_ghost_nodes[jpart];
      ofs << std::setw(idt) << nb_total_edges[jpart];
      ofs << std::setw(idt) << nb_owned_edges[jpart];
      ofs << std::setw(idt) << nb_ghost_edges[jpart];
      ofs << "\n";
    }
    ofs.close();
  }
}

// ------------------------------------------------------------------

// C wrapper interfaces to C++ routines
void atlas__write_load_balance_report (Mesh* mesh, char* filename)
{
  write_load_balance_report( *mesh, std::string(filename) );
}

// ------------------------------------------------------------------


} // namespace actions
} // namespace atlas
