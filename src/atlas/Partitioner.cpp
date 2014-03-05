#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cmath>
#include "Partitioner.hpp"
#include "Mesh.hpp"
#include "FunctionSpace.hpp"
#include "Field.hpp"
#include "Parameters.hpp"
namespace atlas {

Partitioner::Partitioner()
{
}

void Partitioner::partition(Mesh& mesh, int nb_partitions)
{
  FunctionSpace& nodes_2d = mesh.function_space("nodes_2d");
  FieldT<int>& nodes_proc = nodes_2d.field<int>("proc");

  int nb_nodes = nodes_proc.bounds()[1];
  int nodes_per_partition = std::ceil( static_cast<double>(nb_nodes)/static_cast<double>(nb_partitions));
  int partition=0;
  for (int n=0; n<nodes_proc.bounds()[1]; ++n)
  {
    if (n < partition*nodes_per_partition)
    {
      nodes_proc(n) = partition-1;
    }
    else
    {
      ++partition;
      nodes_proc(n) = partition-1;
    }
  }

  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    FunctionSpace& elements = mesh.function_space(f);
    if (elements.metadata<int>("type") == ELEMS)
    {
      FieldT<int>& elem_nodes = elements.field<int>("nodes");
      FieldT<int>& elem_proc = elements.field<int>("proc");
      int nb_elems = elements.bounds()[1];
      int nb_nodes_per_elem = elem_nodes.bounds()[0];
      for (int elem=0; elem<nb_elems; ++elem)
      {
        int proc_min = -1;
        for (int n=0; n<nb_nodes_per_elem; ++n)
        {
          int node = elem_nodes(n,elem);
          if (proc_min<0) proc_min = nodes_proc(node);
          else proc_min = std::min(proc_min, nodes_proc(node) );
        }
        elem_proc(elem) = proc_min;
      }
    }
  }
}

/////////////////////


Partitioner* atlas__Partitioner__new () {
  return new Partitioner();
}

void atlas__Partitioner__delete (Partitioner* This) {
  delete This;
}

/////////////////////

}
