/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iomanip>
#include <fstream>
#include "eckit/filesystem/PathName.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/mpi/mpi.h"
#include "atlas/Mesh.h"
#include "atlas/Nodes.h"
#include "atlas/actions/WriteLoadBalanceReport.h"
#include "atlas/util/IsGhost.h"
#include "atlas/util/IndexView.h"

using atlas::util::IsGhost;

namespace atlas {
namespace actions {

void write_load_balance_report( const Mesh& mesh, const std::string& filename )
{
  std::ofstream ofs;
  if( eckit::mpi::rank() == 0 )
  {
    eckit::PathName path(filename);
    int idt = 10;
    ofs.open( path.localPath(), std::ofstream::out );
  }
  
  write_load_balance_report( mesh, ofs );

  if( eckit::mpi::rank() == 0 )
  {
    ofs.close();
  }
}


void write_load_balance_report( const Mesh& mesh, std::ostream& ofs )
{
  int npart = eckit::mpi::size();
  int root = 0;

  std::vector<int> nb_total_nodes(npart,0);
  std::vector<int> nb_owned_nodes(npart,0);
  std::vector<int> nb_ghost_nodes(npart,0);
  std::vector<double> ghost_ratio_nodes(npart,0);

  std::vector<int> nb_total_edges(npart,0);
  std::vector<int> nb_owned_edges(npart,0);
  std::vector<int> nb_ghost_edges(npart,0);
  std::vector<double> nb_ghost_ratio_edges(npart,0);

  {
    const Nodes& nodes = mesh.nodes();
    IsGhost is_ghost(nodes);
    int nb_nodes = nodes.size();
    int nowned(0);
    int nghost(0);
    for( int n=0; n<nb_nodes; ++n )
    {
      if( is_ghost(n) )
        ++nghost;
      else
        ++nowned;
    }
    ECKIT_MPI_CHECK_RESULT( MPI_Gather( &nb_nodes, 1, MPI_INT,
                                  nb_total_nodes.data(), 1, MPI_INT,
                                  root, eckit::mpi::comm() ) );
    ECKIT_MPI_CHECK_RESULT( MPI_Gather( &nowned, 1, MPI_INT,
                                  nb_owned_nodes.data(), 1, MPI_INT,
                                  root, eckit::mpi::comm() ) );
    ECKIT_MPI_CHECK_RESULT( MPI_Gather( &nghost, 1, MPI_INT,
                                  nb_ghost_nodes.data(), 1, MPI_INT,
                                  root, eckit::mpi::comm() ) );
                                  
    for( size_t p=0; p<npart; ++p )
    {
      ghost_ratio_nodes[p] = static_cast<double>(nb_ghost_nodes[p])/static_cast<double>(nb_owned_nodes[p]);
    }
  }

  bool has_edges = mesh.has_function_space("edges");

  if( has_edges )
  {
    const Nodes& nodes = mesh.nodes();
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
    ECKIT_MPI_CHECK_RESULT( MPI_Gather( &nb_edges, 1, MPI_INT,
                                  nb_total_edges.data(), 1, MPI_INT,
                                  root, eckit::mpi::comm() ) );
    ECKIT_MPI_CHECK_RESULT( MPI_Gather( &nowned, 1, MPI_INT,
                                  nb_owned_edges.data(), 1, MPI_INT,
                                  root, eckit::mpi::comm() ) );
    ECKIT_MPI_CHECK_RESULT( MPI_Gather( &nghost, 1, MPI_INT,
                                  nb_ghost_edges.data(), 1, MPI_INT,
                                  root, eckit::mpi::comm() ) );
  }

  if( eckit::mpi::rank() == 0 )
  {
    int idt = 10;
    ofs << "# STATISTICS\n";
    ofs << std::setw(1)  << "#" << std::setw(5) << "";
    ofs << std::setw(idt) << "nodes";
    ofs << std::setw(idt) << "owned";
    ofs << std::setw(idt) << "ghost";
    ofs << std::setw(idt) << "ratio(%)";
    if( has_edges )
    {
    ofs << std::setw(idt) << "edges";
    ofs << std::setw(idt) << "oedges";
    ofs << std::setw(idt) << "gedges";
    }
    ofs << "\n";
    ofs << std::setw(6)  << "# tot ";
    ofs << std::setw(idt) << std::accumulate(nb_total_nodes.data(),nb_total_nodes.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_owned_nodes.data(),nb_owned_nodes.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_ghost_nodes.data(),nb_ghost_nodes.data()+npart,0);
    ofs << std::setw(idt) << "/";
    if( has_edges )
    {
    ofs << std::setw(idt) << std::accumulate(nb_total_edges.data(),nb_total_edges.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_owned_edges.data(),nb_owned_edges.data()+npart,0);
    ofs << std::setw(idt) << std::accumulate(nb_ghost_edges.data(),nb_ghost_edges.data()+npart,0);
    }
    ofs << "\n";
    ofs << std::setw(6)  << "# max ";
    ofs << std::setw(idt) << *std::max_element(nb_total_nodes.data(),nb_total_nodes.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_owned_nodes.data(),nb_owned_nodes.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_ghost_nodes.data(),nb_ghost_nodes.data()+npart);
    ofs << std::setw(idt) << std::setw(idt) << std::fixed << std::setprecision(2) << *std::max_element(ghost_ratio_nodes.data(),ghost_ratio_nodes.data()+npart) * 100.;
    if( has_edges )
    {
    ofs << std::setw(idt) << *std::max_element(nb_total_edges.data(),nb_total_edges.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_owned_edges.data(),nb_owned_edges.data()+npart);
    ofs << std::setw(idt) << *std::max_element(nb_ghost_edges.data(),nb_ghost_edges.data()+npart);
    }
    ofs << "\n";
    ofs << std::setw(6)  << "# min ";
    ofs << std::setw(idt) << *std::min_element(nb_total_nodes.data(),nb_total_nodes.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_owned_nodes.data(),nb_owned_nodes.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_ghost_nodes.data(),nb_ghost_nodes.data()+npart);
    ofs << std::setw(idt) << std::fixed << std::setprecision(2) << *std::min_element(ghost_ratio_nodes.data(),ghost_ratio_nodes.data()+npart) * 100.;
    if( has_edges )
    {
    ofs << std::setw(idt) << *std::min_element(nb_total_edges.data(),nb_total_edges.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_owned_edges.data(),nb_owned_edges.data()+npart);
    ofs << std::setw(idt) << *std::min_element(nb_ghost_edges.data(),nb_ghost_edges.data()+npart);
    }
    ofs << "\n";
    ofs << std::setw(6)  << "# avg ";
    ofs << std::setw(idt) << std::accumulate(nb_total_nodes.data(),nb_total_nodes.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_owned_nodes.data(),nb_owned_nodes.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_ghost_nodes.data(),nb_ghost_nodes.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::fixed << std::setprecision(2) << std::accumulate(ghost_ratio_nodes.data(),ghost_ratio_nodes.data()+npart,0.)/static_cast<double>(npart) *100.;
    if( has_edges )
    {
    ofs << std::setw(idt) << std::accumulate(nb_total_edges.data(),nb_total_edges.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_owned_edges.data(),nb_owned_edges.data()+npart,0)/npart;
    ofs << std::setw(idt) << std::accumulate(nb_ghost_edges.data(),nb_ghost_edges.data()+npart,0)/npart;
    }
    ofs << "\n";
    ofs << "#----------------------------------------------------\n";
    ofs << "# PER TASK\n";
    ofs << std::setw(6)  << "# part";
    ofs << std::setw(idt) << "nodes";
    ofs << std::setw(idt) << "owned";
    ofs << std::setw(idt) << "ghost";
    ofs << std::setw(idt) << "ratio(%)";
    if( has_edges )
    {
    ofs << std::setw(idt) << "edges";
    ofs << std::setw(idt) << "oedges";
    ofs << std::setw(idt) << "gedges";
    }
    ofs << "\n";
    for( int jpart=0; jpart<npart; ++jpart )
    {
      ofs << std::setw(6)  << jpart;
      ofs << std::setw(idt) << nb_total_nodes[jpart];
      ofs << std::setw(idt) << nb_owned_nodes[jpart];
      ofs << std::setw(idt) << nb_ghost_nodes[jpart];
      ofs << std::setw(idt) << std::fixed << std::setprecision(2) << ghost_ratio_nodes[jpart]*100.;
      if( has_edges )
      {
      ofs << std::setw(idt) << nb_total_edges[jpart];
      ofs << std::setw(idt) << nb_owned_edges[jpart];
      ofs << std::setw(idt) << nb_ghost_edges[jpart];
      }
      ofs << "\n";
    }
  }
}

// ------------------------------------------------------------------

// C wrapper interfaces to C++ routines
void atlas__write_load_balance_report (Mesh* mesh, char* filename)
{
  ATLAS_ERROR_HANDLING( write_load_balance_report( *mesh, std::string(filename) ) );
}

// ------------------------------------------------------------------


} // namespace actions
} // namespace atlas
