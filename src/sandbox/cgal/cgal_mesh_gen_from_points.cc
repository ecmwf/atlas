#include <fstream>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>
#include <CGAL/IO/Polyhedron_VRML_2_ostream.h>

#include "atlas/Parameters.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/Field.hpp"
#include "atlas/FunctionSpace.hpp"

#include "atlas/Gmsh.hpp"

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
typedef GT::Vector_3 Vector_3;
typedef CGAL::Polyhedron_3<GT> Polyhedron;
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
  const double radius_bound   = 0.05;
  const double distance_bound = 0.05;

  CGAL::Surface_mesh_default_criteria_3<Tr> criteria( angular_bound, radius_bound, distance_bound );

  // meshing sphere surface

//  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  /// @todo try out CGAL::Manifold_tag()
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());

  // Output the 2D complex to an OFF file.

//  std::ofstream out("out.off");
//  CGAL::output_surface_facets_to_off (out, c2t3);

  // fill-in Atlas mesh structure

  Mesh* mesh = new Mesh();

  /* nodes */

  const size_t nb_nodes = tr.number_of_vertices();
  std::cout << "nb_nodes = " << nb_nodes << std::endl;

  std::vector<int> bounds(2);
  bounds[0] = Field::UNDEF_VARS;
  bounds[1] = nb_nodes;

  FunctionSpace& nodes    = mesh->add_function_space( new FunctionSpace( "nodes", "Lagrange_P0", bounds ) );

  nodes.metadata().set("type",static_cast<int>(Entity::NODES));

  FieldT<double>& coords = nodes.create_field<double>("coordinates",3);

  std::map< Tr::Vertex_handle, int> vidx;
  int inode = 0;

  for( Tr::Finite_vertices_iterator v = tr.finite_vertices_begin(); v != tr.finite_vertices_end(); ++v)
  {
      vidx[v] = inode;

      const Tr::Point& p = v->point();

      coords(XX,inode) = p.x();
      coords(YY,inode) = p.y();
      coords(ZZ,inode) = p.z();

      ++inode;
//      const double r = std::sqrt( p.x()*p.x() + p.y()*p.y() + p.z()*p.z() );
//      std::cout << vidx[v] << " " <<  p  << std::endl;
  }

  assert( inode == nb_nodes );

  /* triangles */

  const size_t nb_triags = CGAL::Surface_mesher::number_of_facets_on_surface(tr);
  std::cout << "nb_triags = " << nb_triags << std::endl;

  bounds[1] = nb_triags;

  FunctionSpace& triags  = mesh->add_function_space( new FunctionSpace( "triags", "Lagrange_P1", bounds ) );
  triags.metadata().set("type",static_cast<int>(Entity::ELEMS));

  FieldT<int>& triag_nodes   = triags.create_field<int>("nodes",3);

  Point_3 origin (CGAL::ORIGIN);

  size_t tidx = 0;
  for( Tr::Finite_facets_iterator f = tr.finite_facets_begin(); f != tr.finite_facets_end(); ++f )
  {
    const Tr::Cell_handle cell = f->first;
    const int& index = f->second;
    if( cell->is_facet_on_surface(index) == true )
    {
        Tr::Vertex_handle v0 = cell->vertex(tr.vertex_triple_index(index, 0));
        Tr::Vertex_handle v1 = cell->vertex(tr.vertex_triple_index(index, 1));
        Tr::Vertex_handle v2 = cell->vertex(tr.vertex_triple_index(index, 2));

        int idx0 = vidx[ v0 ];
        int idx1 = vidx[ v1 ];
        int idx2 = vidx[ v2 ];

        /* ensure outward pointing normal */

        Vector_3 p0 ( origin, v0->point() );
        Vector_3 n  = CGAL::normal( v0->point(), v1->point(), v2->point() );

        FT innerp = n * p0;

        if( innerp < 0 ) // need to swap an edge of the triag
            std::swap( idx1, idx2 );

        /* define the triag */

        triag_nodes(0,tidx) = F_IDX(idx0);
        triag_nodes(1,tidx) = F_IDX(idx1);
        triag_nodes(2,tidx) = F_IDX(idx2);

        ++tidx;

        //        std::cout << idx0 << " " << idx1 << " " << idx2  << std::endl;
    }
  }

  assert( tidx == nb_triags );

  atlas::Gmsh::write3dsurf(*mesh, std::string("earth.msh") );

  // save in off format
#if 0
   std::ofstream fout("earth.off");
   CGAL::output_surface_facets_to_off( fout, c2t3 );
#endif

   // transform to polyhedron & save to vrml
#if 0
   Polyhedron P;
   CGAL::output_surface_facets_to_polyhedron(c2t3, P);
   std::cin >> P;
   std::ofstream wout("earth.wrl");
//   CGAL::VRML_2_ostream vrml_out( wout );
   CGAL::VRML_1_ostream vrml_out( wout );
   vrml_out << P;
#endif


}
