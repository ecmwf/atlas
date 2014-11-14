/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Util_h
#define atlas_Util_h

#include <cmath>
#include <sstream>

#include "atlas/mpl/MPL.h"
#include "atlas/atlas_config.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Debug.h"
#include "atlas/Parameters.h"

namespace atlas {

inline int microdeg( const double& deg )
{
	return static_cast<int>(deg*1.e6);
}

enum AngleUnit{ DEG=0, RAD=1 };

void colat_to_lat_hemisphere(const int N, const double colat[], double lats[], const AngleUnit unit);

void predict_gaussian_colatitudes_hemisphere(const int N, double colat[]);

void predict_gaussian_latitudes_hemisphere(const int N, double lat[]);


class Flags
{
public:
	static void reset(int& flags, int bit = 0)
	{
		flags = bit;
	}

	static void set(int& flags, int bit)
	{
		flags |= bit;
	}

	static void unset(int& flags, int bit)
	{
		flags &= (~bit);
	}

	static void toggle(int& flags, int bit)
	{
		flags ^= bit;
	}

	static bool check(int flags, int bits)
	{
		return (flags & bits) == bits;
	}

	static bool check_all(int flags, int bits)
	{
		return (flags & bits) == bits;
	}

	static bool check_any(int flags, int bits)
	{
		return flags & bits;
	}

	static std::string bitstr(int flags)
	{
	  char str[9] = {0};
	  int i;
	  for (i=7; i>=0; i--) {
	    str[i] = (flags&1)?'1':'0';
	    flags >>= 1;
	  }
		return std::string(str,9);
	}
};

class Topology : public Flags
{
public:
  enum {
    NONE     = 0,
    GHOST    = (1<<1),
    PERIODIC = (1<<2),
    BC       = (1<<3),
    WEST     = (1<<4),
    EAST     = (1<<5),
    NORTH    = (1<<6),
    SOUTH    = (1<<7)
  };
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
		int i1 = (y/100+NORTH/100*2) >>9;
		int i2 = (x/100+EAST/100*5)  >>10;
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
private:
	static int WEST;
	static int EAST;
	static int NORTH;
	static int SOUTH;
};

class PeriodicTransform
{
private:
	double x_translation_;

public:
	PeriodicTransform()
	{
		x_translation_ = 360.;
	}

	void operator()(double source[2], double dest[2], double direction, double scale = 1.) const
	{
		dest[0] = source[0] + direction*x_translation_*scale;
		dest[1] = source[1];
	}

	void operator()(int source[2], int dest[2], int direction, int scale = 1) const
	{
		dest[0] = source[0] + direction*static_cast<int>(x_translation_*scale);
		dest[1] = source[1];
	}

	void operator()(double inplace[2], double direction, double scale = 1.) const
	{
		inplace[0] = inplace[0] + direction*x_translation_*scale;
		inplace[1] = inplace[1];
	}

	void operator()(int inplace[2], int direction, int scale = 1) const
	{
		inplace[0] = inplace[0] + direction*static_cast<int>(x_translation_*scale);
		inplace[1] = inplace[1];
	}

	void operator()(LatLonPoint& inplace, int direction) const
	{
		inplace.x = inplace.x + direction*static_cast<int>(x_translation_*1.e6);
		inplace.y = inplace.y;
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
		int nb_elem_nodes = elem_nodes.shape(0);
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
