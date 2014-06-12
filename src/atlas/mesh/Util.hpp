/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/mesh/Parameters.hpp"

#include <unistd.h>

#define LOG_DEBUG(WHAT,RANK) \
  if( RANK == -1 || MPL::rank() == RANK) \
    std::cout << "["<< MPL::rank() << "] " << WHAT << std::endl;

#define LOG_DEBUG_VAR(VAR,RANK) \
  if( RANK < 0 || MPL::rank() == RANK) \
    std::cout << "["<< MPL::rank() << "] " << #VAR << " = " << VAR << std::endl;

#define PLOG_DEBUG(WHAT) \
  MPI_Barrier(MPI_COMM_WORLD);\
  LOG_DEBUG(WHAT,-1);\
  MPI_Barrier(MPI_COMM_WORLD); sleep(1);

#define PLOG_DEBUG_VAR(VAR) \
  MPI_Barrier(MPI_COMM_WORLD);\
  LOG_DEBUG_VAR(VAR,-1);\
  MPI_Barrier(MPI_COMM_WORLD); sleep(1);


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
    part    = ArrayView<int,1> (nodes.field("partition") );
    loc_idx = IndexView<int,1> (nodes.field("remote_idx") );
    mypart  = MPL::rank();
  }
  bool operator()(int idx)
  {
    if( part   [idx] != mypart ) return true;
    if( loc_idx[idx] != idx    ) return true;
    return false;
  }
  int mypart;
  ArrayView<int,1> part;
  IndexView<int,1> loc_idx;
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
