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

  ArrayView<int,1> nloen       = t_->nloen();
  ArrayView<int,1> n_regions   = t_->n_regions();
  ArrayView<int,1> nfrstlat    = t_->nfrstlat();
  ArrayView<int,1> nlstlat     = t_->nlstlat();
  ArrayView<int,1> nptrfrstlat = t_->nptrfrstlat();
  ArrayView<int,2> nsta        = t_->nsta();
  ArrayView<int,2> nonl        = t_->nonl();


  int i(0);
  int maxind(0);
  std::vector<int> iglobal(t_->ndgl()*nlonmax,-1);

  for( int jgl=0; jgl<t_->ndgl(); ++jgl )
  {
    for( int jl=0; jl<nloen[jgl]; ++jl )
    {
      ++i;
      iglobal[jgl*nlonmax+jl] = i;
      maxind = std::max(maxind,jgl*nlonmax+jl);
    }
  }
  int iproc(0);
  for( int ja=0; ja<t_->n_regions_NS(); ++ja )
  {
    for( int jb=0; jb<n_regions[ja]; ++jb )
    {
      for( int jgl=nfrstlat[ja]-1; jgl<nlstlat[ja]; ++jgl )
      {
        int igl = nptrfrstlat[ja] + jgl - nfrstlat[ja];
        for( int jl=nsta(jb,igl)-1; jl<nsta(jb,igl)+nonl(jb,igl)-1; ++jl )
        {
          int ind = iglobal[jgl*nlonmax+jl] - 1;
          ASSERT( ind >= 0 );
          if( ind >= grid().npts() ) throw eckit::OutOfRange(ind,grid().npts(),Here());
          part[ind] = iproc;
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

