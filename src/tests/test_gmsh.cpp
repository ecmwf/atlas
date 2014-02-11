#include "datastruct/Gmsh.hpp"
#include "datastruct/BuildEdges.hpp"
#include "datastruct/BuildDualMesh.hpp"
#include "datastruct/BuildPeriodicBoundaries.hpp"
#include "datastruct/Partitioner.hpp"
#include <mpi.h>

using namespace ecmwf;
int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);

  Mesh& mesh = Gmsh::read("T47.msh");
  //Mesh& mesh = Gmsh::read("untitled.msh");
  //Mesh& mesh = Gmsh::read("test_no_edges.msh");

  //Partitioner::partition(mesh,1);

  build_periodic_boundaries(mesh);
  build_edges(mesh);
  build_dual_mesh(mesh);
  Gmsh::write(mesh,"bla.msh");
  
  MPI_Finalize();
  return 0;
}
