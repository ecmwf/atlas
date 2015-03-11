/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/mpi/mpi.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/grids/ReducedGrid.h"
#include "atlas/trans/TransPartitioner.h"
#include "atlas/trans/Trans.h"

namespace atlas {
namespace trans {

TransPartitioner::TransPartitioner( const grids::ReducedGrid& grid, const Trans& trans ) :
  Partitioner(grid),
  t_( const_cast<Trans*>(&trans)), owned_(false)
{
  ASSERT( t_ != NULL );
  Partitioner::set_nb_partition( t_->nproc() );
}

TransPartitioner::TransPartitioner( const grids::ReducedGrid& grid ) :
  Partitioner(grid),
  t_( new Trans(grid,0) ), owned_(true)
{
  ASSERT( t_ != NULL );
  Partitioner::set_nb_partition( t_->nproc() );
}

TransPartitioner::~TransPartitioner()
{
  if( owned_ )
    delete t_;
}

void TransPartitioner::partition(int part[]) const
{
  if( dynamic_cast<const grids::ReducedGrid*>(&grid()) == NULL )
    throw eckit::BadCast("Grid is not a grids::ReducedGrid type. Cannot partition using IFS trans",Here());

  int nlonmax = dynamic_cast<const grids::ReducedGrid*>(&grid())->nlonmax();

  int i(0);
  std::vector<int> iglobal(t_->ndgl()*nlonmax);
  for( int jgl=1; jgl<=t_->ndgl(); ++jgl )
  {
    for( int jl=1; jl<=t_->nloen()[jgl-1]; ++jl )
    {
      ++i;
      iglobal[(jgl-1)*nlonmax+(jl-1)] = i;
    }
  }

  const int stride=t_->ndgl()+t_->n_regions_NS()-1;

  int iproc(0);
  for( int ja=1; ja<=t_->n_regions_NS(); ++ja )
  {
    for( int jb=1; jb<=t_->n_regions()[ja-1]; ++jb )
    {
      for( int jgl=t_->nfrstlat()[ja-1]; jgl<=t_->nlstlat()[ja-1]; ++jgl )
      {
        int igl = t_->nptrfrstlat()[ja-1] + jgl - t_->nfrstlat()[ja-1];

        int idx = (jb-1)*stride+(igl-1);
        for( int jl=t_->nsta()[idx]; jl<=t_->nsta()[idx]+t_->nonl()[idx]-1; ++jl )
        {
          int ind = iglobal[(jgl-1)*nlonmax+(jl-1)];
          part[ind-1] = iproc;
        }
      }
      ++iproc;
    }
  }
}

int TransPartitioner::nb_bands() const
{
  ASSERT( t_!= NULL );
  return t_->n_regions_NS();
}

int TransPartitioner::nb_regions(int b) const
{
  ASSERT( t_!= NULL );
  return t_->n_regions()[b];
}

}
}

