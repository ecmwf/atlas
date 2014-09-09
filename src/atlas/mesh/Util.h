/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_mesh_Util_h
#define atlas_mesh_Util_h

#include <cmath>
#include <sstream>

#include "atlas/mpl/MPL.h"
#include "atlas/atlas_config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/FunctionSpace.h"
#include "atlas/mesh/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Debug.h"
#include "atlas/mesh/Parameters.h"

namespace atlas {

inline int microdeg( const double& deg )
{
	return static_cast<int>(deg*1.e6);
}

struct BC
{
	static int WEST;
	static int EAST;
	static int NORTH;
	static int SOUTH;
};

struct LatLonPoint
{
	// Storage is in microdegrees
	// This structure is used in sorting algorithms, and uses less memory than
	// if x and y were in double precision.
	LatLonPoint() {}
	LatLonPoint( int x_, int y_ )
	{
		x = x_;
		y = y_;
	}
	LatLonPoint( double x_, double y_ )
	{
		x = microdeg(x_);
		y = microdeg(y_);
	}
	LatLonPoint( const double coord[2] )
	{
		x = microdeg(coord[XX]);
		y = microdeg(coord[YY]);
	}
	LatLonPoint( const ArrayView<int,1>& coord )
	{
		x = coord[XX];
		y = coord[YY];
	}
	LatLonPoint( const ArrayView<double,1>& coord )
	{
		x = microdeg(coord[XX]);
		y = microdeg(coord[YY]);
	}

	int uid() const
	{
		int i1 = (y+BC::NORTH*2) >>9;
		int i2 = (x+BC::EAST*5)  >>10;
		ASSERT( i1 > 0);
		ASSERT( i2 > 0);
		int pow = 10;
		while(i2 >= pow)
				pow *= 10;
		int id = i1*pow + i2;
		ASSERT( id > 0 );
		return id;
	}

	mutable int x, y;
	bool operator < (const LatLonPoint& other) const
	{
		if( y > other.y  ) return true;
		if( y == other.y ) return (x < other.x);
		return false;
	}
};

struct IsGhost
{
	IsGhost( FunctionSpace& nodes )
	{
		part_	 = ArrayView<int,1> (nodes.field("partition") );
		ridx_	 = IndexView<int,1> (nodes.field("remote_idx") );
		mypart_ = MPL::rank();
	}
	IsGhost( FunctionSpace& nodes, int mypart )
	{
		part_	 = ArrayView<int,1> (nodes.field("partition") );
		ridx_	 = IndexView<int,1> (nodes.field("remote_idx") );
		mypart_ = mypart;
	}

	bool operator()(int idx)
	{
		if( part_[idx] != mypart_ ) return true;
		if( ridx_[idx] != idx     ) return true;
		return false;
	}
	int mypart_;
	ArrayView<int,1> part_;
	IndexView<int,1> ridx_;
};

struct ComputeUniqueElementIndex
{
	ComputeUniqueElementIndex( const FunctionSpace& nodes )
	{
		coords = ArrayView<double,2> ( nodes.field("coordinates") );
	}

	int operator()( const IndexView<int,1>& elem_nodes ) const
	{
		double centroid[2];
		centroid[XX] = 0.;
		centroid[YY] = 0.;
		int nb_elem_nodes = elem_nodes.extents()[0];
		for( int jnode=0; jnode<nb_elem_nodes; ++jnode )
		{
			centroid[XX] += coords( elem_nodes(jnode), XX );
			centroid[YY] += coords( elem_nodes(jnode), YY );
		}
		centroid[XX] /= static_cast<double>(nb_elem_nodes);
		centroid[YY] /= static_cast<double>(nb_elem_nodes);
		return LatLonPoint( centroid[XX], centroid[YY] ).uid();
	}

	ArrayView<double,2> coords;
};


struct Face
{
	ElementRef& operator[](int i) { return elems[i]; }
	bool is_bdry() const { return (elems[1].f < 0); }
	ElementRef elems[2];
};

void accumulate_faces(
		FunctionSpace& func_space,
		std::vector< std::vector<int> >& node_to_face,
		std::vector<int>& face_nodes_data,
		std::vector< Face >& connectivity_edge_to_elem,
		int& nb_faces,
		int& nb_inner_faces );

} // namespace atlas

#endif
