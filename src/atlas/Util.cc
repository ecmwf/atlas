/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Util.h"

namespace atlas {

int LatLonPoint::WEST  = 0;
int LatLonPoint::EAST  = microdeg( 360. );
int LatLonPoint::NORTH = microdeg(  90. );
int LatLonPoint::SOUTH = microdeg( -90. );


template<>
gidx_t LatLonPoint::uidT<int>() const
{
  return uid32();
}

template<>
gidx_t LatLonPoint::uidT<long int>() const
{
  return uid64();
}

gidx_t LatLonPoint::uid() const
{
  return uid32();
//  return uidT<gidx_t>();
}

void accumulate_faces(
		FunctionSpace& func_space,
		std::vector< std::vector<int> >& node_to_face,
		std::vector<int>& face_nodes_data,
		std::vector< Face >& connectivity_edge_to_elem,
		int& nb_faces,
		int& nb_inner_faces )
{
	IndexView<int,2> elem_nodes( func_space.field( "nodes" ) );
	int nb_elems = func_space.shape(0);
	int nb_nodes_in_face = 2;

	std::vector< std::vector<int> > face_node_numbering;
	int nb_faces_in_elem;
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

	for (int e=0; e<nb_elems; ++e)
	{
		for (int f=0; f<nb_faces_in_elem; ++f)
		{
			bool found_face = false;

			std::vector<int> face_nodes(nb_nodes_in_face);
			for (int jnode=0; jnode<nb_nodes_in_face; ++jnode)
				face_nodes[jnode] = elem_nodes(e,face_node_numbering[f][jnode]);

			int node = face_nodes[0];
			for( int jface=0; jface< node_to_face[node].size(); ++jface )
			{
				int face = node_to_face[node][jface];
				int nb_matched_nodes = 0;
				if (nb_nodes_in_face>1) // 2D or 3D
				{
					for( int jnode=0; jnode<nb_nodes_in_face; ++jnode)
					{
						int other_node = face_nodes[jnode];
						for( int iface=0; iface<node_to_face[other_node].size(); ++iface )
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
				for (int n=0; n<nb_nodes_in_face; ++n)
				{
					node_to_face[face_nodes[n]].push_back(nb_faces);
					face_nodes_data.push_back(face_nodes[n]);
				}
				++nb_faces;
			}
		}
	}
}

namespace {
static const double rad_to_deg = 180.*M_1_PI;
}

void colat_to_lat_hemisphere(const int N, const double colat[], double lats[], const AngleUnit unit)
{
  std::copy( colat, colat+N, lats );
  double pole = (unit == DEG ? 90. : M_PI_2);
  for (int i=0; i<N; ++i)
    lats[i]=pole-lats[i];
}

void predict_gaussian_colatitudes_hemisphere(const int N, double colat[])
{
  double z;
  for( int i=0; i<N; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*N+2.);
    colat[i] = ( z+1./(tan(z)*(8.*(2.*N)*(2.*N))) ) * rad_to_deg;
  }
}

void predict_gaussian_latitudes_hemisphere(const int N, double lats[])
{
  std::vector<double> colat(N);
  predict_gaussian_colatitudes_hemisphere(N,colat.data());
  colat_to_lat_hemisphere(N,colat.data(),lats,DEG);
}

} // namespace atlas
