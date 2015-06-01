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

//#include "atlas/mpi/mpi.h"
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
	return static_cast<int>(deg*1.e6 + 0.5);
}

enum AngleUnit{ DEG=0, RAD=1 };

void colat_to_lat_hemisphere(const int N, const double colat[], double lats[], const AngleUnit unit);

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


struct LonLatPoint
{
	// Storage is in microdegrees
	// This structure is used in sorting algorithms, and uses less memory than
	// if x and y were in double precision.
	LonLatPoint() {}
	LonLatPoint( int x_, int y_ )
	{
		x = x_;
		y = y_;
	}
  LonLatPoint( long x_, long y_ )
	{
		x = x_;
		y = y_;
	}
  LonLatPoint( double x_, double y_ )
	{
		x = microdeg(x_);
		y = microdeg(y_);
	}
	LonLatPoint( const double coord[2] )
	{
		x = microdeg(coord[LON]);
		y = microdeg(coord[LAT]);
	}
	LonLatPoint( const ArrayView<int,1>& coord )
	{
		x = coord[LON];
		y = coord[LAT];
	}
	LonLatPoint( const ArrayView<double,1>& coord )
	{
		x = microdeg(coord[LON]);
		y = microdeg(coord[LAT]);
	}

  long uid64() const;

  int uid32() const;

  gidx_t uid() const;

	mutable int x, y;
	bool operator < (const LonLatPoint& other) const
	{
		if( y > other.y  ) return true;
		if( y == other.y ) return (x < other.x);
		return false;
	}
private:

  template<typename T> gidx_t uidT() const;

  static int WEST;
	static int EAST;
	static int NORTH;
	static int SOUTH;
};

template<> inline gidx_t LonLatPoint::uidT<int >() const { return uid32(); }
template<> inline gidx_t LonLatPoint::uidT<long>() const { return uid64(); }

inline int LonLatPoint::uid32() const
{
  // max precision is 0.02 degree
  int iy = static_cast<int>((2*NORTH-y)*5e-5);
  int ix = static_cast<int>((x+2*EAST)*5e-5);
  iy <<= 17;
  int id = iy | ix;
  return id;
}

inline long LonLatPoint::uid64() const
{
  // max precision is 1 microdegree
  long iy = static_cast<long>((4*NORTH-y));
  long ix = static_cast<long>((x+4*EAST));
  iy <<= 31;
  long id = iy | ix;
  return id;
}

inline gidx_t LonLatPoint::uid() const {
  return uidT<gidx_t>();
}


class PeriodicTransform
{
protected:
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

	void operator()(LonLatPoint& inplace, int direction) const
	{
		inplace.x = inplace.x + direction*microdeg(x_translation_);
		inplace.y = inplace.y;
	}

};


struct IsGhost
{
	IsGhost( FunctionSpace& nodes )
	{
		part_	 = ArrayView<int,1> (nodes.field("partition") );
		ridx_	 = IndexView<int,1> (nodes.field("remote_idx") );
		mypart_ = eckit::mpi::rank();
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


struct ComputeUid
{
  ComputeUid() {}
	ComputeUid( const FunctionSpace& nodes ):
    funcspace(&nodes)
	{
		update();
	}

  gidx_t operator()( const double crd[] ) const
	{
		return LonLatPoint( crd ).uid();
	}

	gidx_t operator()( int node ) const
	{
		return LonLatPoint( lonlat[node] ).uid();
	}

  gidx_t operator()( const IndexView<int,1>& elem_nodes ) const
	{
		double centroid[2];
		centroid[LON] = 0.;
		centroid[LAT] = 0.;
		int nb_elem_nodes = elem_nodes.shape(0);
		for( int jnode=0; jnode<nb_elem_nodes; ++jnode )
		{
			centroid[LON] += lonlat( elem_nodes(jnode), LON );
			centroid[LAT] += lonlat( elem_nodes(jnode), LAT );
		}
		centroid[LON] /= static_cast<double>(nb_elem_nodes);
		centroid[LAT] /= static_cast<double>(nb_elem_nodes);
		return LonLatPoint( centroid[LON], centroid[LAT] ).uid32();
	}

  gidx_t operator()( double crds[], int npts ) const
	{
		double centroid[2];
		centroid[LON] = 0.;
		centroid[LAT] = 0.;
		for( int jnode=0; jnode<npts; ++jnode )
		{
			centroid[LON] += crds[jnode*2+LON];
			centroid[LAT] += crds[jnode*2+LAT];
		}
		centroid[LON] /= static_cast<double>(npts);
		centroid[LAT] /= static_cast<double>(npts);
		return LonLatPoint( centroid[LON], centroid[LAT] ).uid32();
	}


  void update()
  {
    lonlat = ArrayView<double,2> ( funcspace->field("lonlat") );
  }
private:
  const FunctionSpace* funcspace;
  ArrayView<double,2> lonlat;
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
