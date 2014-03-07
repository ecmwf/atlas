#include <fstream>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include "atlas/Parameters.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/Field.hpp"
#include "atlas/FunctionSpace.hpp"

// #define ER   6371     // earth readius
#define ER   1     // earth readius
#define ER2  ER*ER    // earth readius squared
#define DER2 2*ER2    // bouding sphere squared radius

// default triangulation for Surface_mesher

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

FT sphere_function (Point_3 p)
{
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2 + y2 + z2 - ER2 ;
}

using namespace atlas;

int main()
{
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface

  Sphere_3 bound_sphere(CGAL::ORIGIN, DER2 );

  Surface_3 surface( sphere_function, bound_sphere);

  // defining meshing criteria

  const double angular_bound  = 30.;
  const double radius_bound   = 0.25;
  const double distance_bound = 0.25;

  CGAL::Surface_mesh_default_criteria_3<Tr> criteria( angular_bound, radius_bound, distance_bound );

  // meshing sphere surface

  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  /// @todo try out CGAL::Manifold_tag()
  //  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());

  // Output the 2D complex to an OFF file.

//  std::ofstream out("out.off");
//  CGAL::output_surface_facets_to_off (out, c2t3);

  // list the vertices

  CGAL::set_ascii_mode( std::cout);

  std::map< Tr::Vertex_handle, int> vidx;
  int inum = 0;

  for( Tr::Finite_vertices_iterator v = tr.finite_vertices_begin(); v != tr.finite_vertices_end(); ++v)
  {
      vidx[v] = inum++;

      const typename Tr::Point& p = v->point();

//      const double r = std::sqrt( p.x()*p.x() + p.y()*p.y() + p.z()*p.z() );

      std::cout << vidx[v] << " " <<  p  << std::endl;
  }

  for( Tr::Finite_facets_iterator f = tr.finite_facets_begin(); f != tr.finite_facets_end(); ++f )
  {
    const typename Tr::Cell_handle cell = f->first;
    const int& index = f->second;
    if( cell->is_facet_on_surface(index) == true )
    {
        const int idx0 = vidx[ cell->vertex(tr.vertex_triple_index(index, 0)) ];
        const int idx1 = vidx[ cell->vertex(tr.vertex_triple_index(index, 1)) ];
        const int idx2 = vidx[ cell->vertex(tr.vertex_triple_index(index, 2)) ];

//        std::cout << idx0 << " " << idx1 << " " << idx2  << std::endl;
    }
  }

  // fill-in Atlas mesh structure

  Mesh* mesh = new Mesh();

  /* nodes */

  const size_t nb_nodes = tr.number_of_vertices();

  std::vector<int> bounds(2);
  bounds[0] = Field::UNDEF_VARS;
  bounds[1] = nb_nodes;

  std::cout << "nb_nodes = " << nb_nodes << std::endl;

  FunctionSpace& nodes_3d    = mesh->add_function_space( new FunctionSpace( "nodes_3d", "Lagrange_P0", bounds ) );

  nodes_3d.metadata().set("type",static_cast<int>(Entity::NODES));
}
