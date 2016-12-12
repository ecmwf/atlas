#include <cmath>
#include <ctime>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm>
#include "eckit/geometry/Point2.h"
#include "atlas/grid/regular/Regular.h"
#include "atlas/grid/partitioners/CheckerBoardPartitioner.h"
#include "atlas/internals/Functions.h"
#include "atlas/parallel/mpi/mpi.h"

using atlas::internals::microdeg;

namespace atlas {
namespace grid {
namespace partitioners {

CheckerBoardPartitioner::CheckerBoardPartitioner(const grid::Grid& grid) :
	Partitioner(grid,eckit::mpi::size())
{
	// defaults
	configure_defaults(grid);
	
	// perform check
	if (checkerboard_ && nparts_%nbands_!=0) throw eckit::BadValue("number of bands doesn't divide number of partitions",Here());

}

CheckerBoardPartitioner::CheckerBoardPartitioner(const grid::Grid& grid, int N) :
  Partitioner(grid,N)
{
	// defaults
	configure_defaults(grid);
	
	// arguments
	nparts_=N;

	// perform check
	if (checkerboard_ && nparts_%nbands_!=0) throw eckit::BadValue("number of bands doesn't divide number of partitions",Here());
}

CheckerBoardPartitioner::CheckerBoardPartitioner(const grid::Grid& grid, int N, int nbands) :
  Partitioner(grid,N)
{
	// defaults
	configure_defaults(grid);
	
	// arguments
	nparts_=N;
	nbands_=nbands;

	// perform check
	if (checkerboard_ && nparts_%nbands_!=0) throw eckit::BadValue("number of bands doesn't divide number of partitions",Here());
}

CheckerBoardPartitioner::CheckerBoardPartitioner(const grid::Grid& grid, int N, int nbands, bool checkerboard) :
  Partitioner(grid,N)
{
	// defaults
	configure_defaults(grid);
	
	// arguments
	nparts_=N;
	nbands_=nbands;
	checkerboard_=checkerboard;

	// perform check
	if (checkerboard_ && nparts_%nbands_!=0) throw eckit::BadValue("number of bands doesn't divide number of partitions",Here());
}

void CheckerBoardPartitioner::configure_defaults(const grid::Grid& grid) {
	// default number of parts
	nparts_=eckit::mpi::size();

	// grid dimensions
  const grid::regular::Regular* rg = dynamic_cast<const grid::regular::Regular*>(&grid);
	nx_=rg->nlonmax();
	ny_=rg->nlat();
	
	// default number of bands
	double zz=sqrt( (double) (nparts_*ny_)/nx_ );		// aim at +/-square regions
	nbands_=(int) zz+0.5;
	if ( nbands_ <1 ) nbands_=1;							// at least one band
	if ( nbands_ > nparts_ ) nbands_=nparts_;	// not more bands than procs
	
	// default checkerboard
	checkerboard_=true;
	
	// true checkerboard means nbands must divide nparts
	if (checkerboard_) {
		while ( nparts_ % nbands_ !=0 ) nbands_--;
	}
}

bool compare_Y_X(const CheckerBoardPartitioner::NodeInt& node1, const CheckerBoardPartitioner::NodeInt& node2)
{
	// comparison of two locations; X1 < X2 if it's to the south, then to the east.
  if( node1.y <  node2.y ) return true;
  if( node1.y == node2.y ) return (node1.x < node2.x);
  return false;
}

bool compare_X_Y(const CheckerBoardPartitioner::NodeInt& node1, const CheckerBoardPartitioner::NodeInt& node2)
{
	// comparison of two locations; X1 < X2 if it's to the east, then to the south.
  if( node1.x <  node2.x ) return true;
  if( node1.x == node2.x ) return (node1.y < node2.y);
  return false;
}

void CheckerBoardPartitioner::partition(int nb_nodes, NodeInt nodes[], int part[]) const
{

  int nparts = nparts_;
	int nbands = nbands_;
	int remainder;
	
  /*
  Sort nodes from south to north (increasing y), and west to east (increasing x). Now we can easily split
  the points in bands. Note this may not be necessary, as it could be
  already by construction in this order, but then sorting is really fast
  */
  
  std::cout << "nbands = " << nbands << std::endl;

  /*
	Number of procs per band
  */
	std::vector< int > npartsb(nbands,0);			// number of procs per band
	remainder=nparts;
	for (int iband=0;iband<nbands;iband++)
	{
		npartsb[iband]=nparts/nbands;
		remainder-=npartsb[iband];
	}
	// distribute remaining procs over first bands
	for (int iband=0;iband<remainder;iband++) ++npartsb[iband];
	
	
  /*
	Number of gridpoints per band
  */
	std::vector< int > ngpb(nbands,0);
	// split latitudes?
	if (true)
	{
		remainder=nb_nodes;
		for (int iband=0;iband<nbands;iband++)
		{
			ngpb[iband]=(nb_nodes*npartsb[iband])/nparts;
			remainder-=ngpb[iband];
		}
		// distribute remaining gridpoints over first bands
		for (int iband=0;iband<remainder;iband++) ++ngpb[iband];

	} else {
		remainder=ny_;
		for (int iband=0;iband<nbands;iband++)
		{
			ngpb[iband]=nx_*(( ny_*npartsb[iband])/nparts);
			remainder-=ngpb[iband]/nx_;
		}
		// distribute remaining rows over first bands
		for (int iband=0;iband<remainder;iband++) ngpb[iband]+=nx_;
	}

	//for (int iband=0;iband<nbands; iband++ ) std::cout << "band " << iband << " : nparts = " << npartsb[iband] << "; ngpb = " << ngpb[iband] << std::endl;
	
	

	// sort nodes according to Y first, to determine bands
  std::sort( nodes, nodes+nb_nodes, compare_Y_X);
  
  //std::cout << __LINE__ << ",  in " << __FUNCTION__ << std::endl;
    
	// for each band, select gridpoints belonging to that band, and sort them according to X first
	int offset=0;
	int jpart=0;
	for (int iband=0;iband<nbands;iband++)
	{
		// sort according to X first
		std::sort( nodes+offset, nodes+offset+ngpb[iband], compare_X_Y);
		
		
		// number of gridpoints per task
		std::vector< int > ngpp(npartsb[iband],0);
		remainder=ngpb[iband];
		//std::cout << "remainder = " << remainder << std::endl;
		for (int ipart=0;ipart<npartsb[iband];ipart++)
		{
			ngpp[ipart]=ngpb[iband]/npartsb[iband];
			remainder-=ngpp[ipart];
			//std::cout << "remainder = " << remainder << std::endl;
		}
		// distribute remaining gridpoints over first parts
		for (int ipart=0;ipart<remainder;ipart++) ++ngpp[ipart];
		
		/*
		std::cout << "ngpb = " << ngpb[iband] << "\n";
		for (int ipart=0;ipart<npartsb[iband];ipart++) std::cout << "part " << ipart << "; ngpp = " << ngpp[ipart] << ". ";
		std::cout << std::endl;
		*/

		// set partition number for each part
		for (int ipart=0;ipart<npartsb[iband];ipart++ )
		{
			//std::cout << "looping over points " << offset << " to " << offset+ngpp[ipart] << std::endl;
			for (int jj=offset;jj<offset+ngpp[ipart];jj++)
			{
				part[nodes[jj].n] = jpart;
			}
			offset+=ngpp[ipart];
			++jpart;
		}

		//if (iband>0) {  sleep(eckit::mpi::size()-eckit::mpi::rank()); ASSERT(0);}

	}	
	
#ifdef gnarls
	// nice output of partition
	if ( eckit::mpi::rank()==0 ) {
		for (int iy=0;iy<ny_;iy++) {
			for (int ix=0;ix<nx_;ix++) {
				std::cout  << std::setw(3) << part[iy*nx_+ix] << " ";
			}
			std::cout << "\n";
		}
		std::cout << std::endl;
	}
  //ASSERT(0);
#endif

}

void CheckerBoardPartitioner::partition(int part[]) const
{
  if( nparts_ == 1 ) // trivial solution, so much faster
  {
    for(size_t j = 0; j < grid().npts(); ++j)
      part[j] = 0;
  }
  else
  {
    std::vector<NodeInt> nodes(grid().npts());
    int n(0);

		for(size_t iy = 0; iy < ny_; ++iy)
		{
			for(size_t ix = 0; ix < nx_; ++ix)
			{
				nodes[n].x = ix;
				nodes[n].y = iy;
				nodes[n].n = n;
				++n;
			}
		}

		partition(grid().npts(), nodes.data(), part);
  }
}

} // namespace partitioners
} // namespace grid
} // namespace atlas

namespace {
    atlas::grid::partitioners::PartitionerBuilder<atlas::grid::partitioners::CheckerBoardPartitioner> __CheckerBoard("CheckerBoard");
}



