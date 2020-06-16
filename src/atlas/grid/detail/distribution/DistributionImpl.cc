/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <ostream>
#include <vector>

#include "DistributionImpl.h"

#include "atlas/grid/Grid.h"
#include "atlas/util/Config.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"

namespace atlas {
namespace grid {

namespace {
std::string distribution_type( int N, const Partitioner& p = Partitioner() ) {
    if ( N == 1 ) {
        return "serial";
    }
    if ( not p ) {
        return "custom";
    }
    return p.type();
}
}  // namespace

DistributionImpl::DistributionImpl( const Grid & grid, const eckit::Parametrisation & config )
{
  bool light = false;

  config.get ("light", light);

  if (light)
    {  
      gridsize_ = grid.size ();
      nb_partitions_ = atlas::mpi::comm ().size ();
      blocksize_ = 1;

      config.get ("blocksize", blocksize_);
     
      nb_blocks_ = gridsize_ / blocksize_;
     
      if (gridsize_ % blocksize_)
        nb_blocks_++;
     
      nb_pts_.reserve (nb_partitions_);
     
      for (idx_t iproc = 0; iproc < nb_partitions_; iproc++)
         {
           // Approximate values
           gidx_t imin = blocksize_ * (((iproc + 0) * nb_blocks_) / nb_partitions_);
           gidx_t imax = blocksize_ * (((iproc + 1) * nb_blocks_) / nb_partitions_);

           while (imin > 0)
             if (partition (imin-blocksize_) == iproc)
               imin -= blocksize_;
             else
               break;

           while (partition (imin) < iproc)
             imin += blocksize_;

           while (partition (imax-1) == iproc + 1)
             imax -= blocksize_;

           while (imax + blocksize_ <= gridsize_)
             if (partition (imax) == iproc)
               imax += blocksize_;
             else
               break;

           imax = std::min (imax, (gidx_t)gridsize_);
           nb_pts_.push_back (imax-imin);
         }
     
      max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
      min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
    }
  else
    {
      Partitioner partitioner (config);
      setupWithPartitioner (grid, partitioner);
    }

}

DistributionImpl::DistributionImpl( const Grid& grid ) :
    nb_partitions_( 1 ),
    part_( grid.size(), 0 ),
    nb_pts_( nb_partitions_, grid.size() ),
    max_pts_( grid.size() ),
    min_pts_( grid.size() ),
    type_( distribution_type( nb_partitions_ ) ) {}

void DistributionImpl::setupWithPartitioner (const Grid & grid, const Partitioner & partitioner)
{
  part_.resize (grid.size ());
  partitioner.partition( grid, part_.data() );
  nb_partitions_ = partitioner.nb_partitions();

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // new
  size_t size     = part_.size();
  int num_threads = atlas_omp_get_max_threads();

  std::vector<std::vector<int> > nb_pts_per_thread( num_threads, std::vector<int>( nb_partitions_ ) );
  atlas_omp_parallel {
      int thread   = atlas_omp_get_thread_num();
      auto& nb_pts = nb_pts_per_thread[thread];
      atlas_omp_for( size_t j = 0; j < size; ++j ) {
          int p = part_[j];
          ++nb_pts[p];
      }
  }

  nb_pts_.resize( nb_partitions_, 0 );
  for ( int thread = 0; thread < num_threads; ++thread ) {
      for ( int p = 0; p < nb_partitions_; ++p ) {
          nb_pts_[p] += nb_pts_per_thread[thread][p];
      }
  }


  // ==============================================
  // previous
  //
  // nb_pts_.resize( nb_partitions_, 0 );
  // for ( idx_t j = 0, size = static_cast<idx_t>( part_.size() ); j < size; ++j ) {
  //     ++nb_pts_[part_[j]];
  // }
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
  min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
  type_    = distribution_type( nb_partitions_, partitioner );
}

DistributionImpl::DistributionImpl( const Grid& grid, const Partitioner& partitioner ) {
  setupWithPartitioner (grid, partitioner);
}

DistributionImpl::DistributionImpl( int nb_partitions, idx_t npts, int part[], int part0 ) {
    part_.assign( part, part + npts );
    if ( nb_partitions == 0 ) {
        std::set<int> partset( part_.begin(), part_.end() );
        nb_partitions_ = static_cast<idx_t>( partset.size() );
    }
    else {
        nb_partitions_ = nb_partitions;
    }
    nb_pts_.resize( nb_partitions_, 0 );
    for ( idx_t j = 0, size = static_cast<idx_t>( part_.size() ); j < size; ++j ) {
        part_[j] -= part0;
        ++nb_pts_[part_[j]];
    }
    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
    type_    = distribution_type( nb_partitions_ );
}

DistributionImpl::DistributionImpl( int nb_partitions, partition_t&& part ) :
    nb_partitions_( nb_partitions ), part_( std::move( part ) ), nb_pts_( nb_partitions_, 0 ) {
    size_t size     = part_.size();
    int num_threads = atlas_omp_get_max_threads();
    std::vector<std::vector<int> > nb_pts_per_thread( num_threads, std::vector<int>( nb_partitions_ ) );
    atlas_omp_parallel {
        int thread   = atlas_omp_get_thread_num();
        auto& nb_pts = nb_pts_per_thread[thread];
        atlas_omp_for( size_t j = 0; j < size; ++j ) {
            int p = part_[j];
            ++nb_pts[p];
        }
    }
    for ( int thread = 0; thread < num_threads; ++thread ) {
        for ( int p = 0; p < nb_partitions_; ++p ) {
            nb_pts_[p] += nb_pts_per_thread[thread][p];
        }
    }

    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
    type_    = distribution_type( nb_partitions_ );
}

DistributionImpl::~DistributionImpl() = default;

void DistributionImpl::print( std::ostream& s ) const {
    s << "Distribution( "
      << "type: " << type_ << ", nb_points: " << part_.size() << ", nb_partitions: " << nb_pts_.size() << ", parts : [";
    for ( idx_t i = 0, size = static_cast<idx_t>( part_.size() ); i < size; i++ ) {
        if ( i != 0 ) {
            s << ',';
        }
        s << part_[i];
    }
    s << ']';
}


DistributionImpl* atlas__GridDistribution__new( idx_t npts, int part[], int part0 ) {
    return new DistributionImpl( 0, npts, part, part0 );
}

void atlas__GridDistribution__delete( DistributionImpl* This ) {
    delete This;
}

void atlas__GridDistribution__nb_pts( DistributionImpl* This, idx_t nb_pts[] ) {
    const auto& nb_pts_ = This->nb_pts();
    std::copy( nb_pts_.begin(), nb_pts_.end(), &nb_pts[0] );
}

idx_t atlas__atlas__GridDistribution__nb_partitions( DistributionImpl* This ) {
    return This->nb_partitions();
}


}  // namespace grid
}  // namespace atlas
