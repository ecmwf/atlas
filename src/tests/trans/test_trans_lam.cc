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
#include <cmath>
#include <iomanip>
#include <numeric>

#include "eckit/config/YAMLConfiguration.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/option/TransOptions.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/local/TransLocal.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"

#include "tests/AtlasTestEnvironment.h"

// #if ATLAS_HAVE_TRANS
// #include "atlas/library/config.h"
// #if ATLAS_HAVE_ECTRANS
// #include "ectrans/transi.h"
// #else
// #include "transi/trans.h"
// #endif
// #endif

using namespace eckit;

using atlas::array::Array;
using atlas::array::ArrayView;
using atlas::array::make_view;

namespace atlas {
namespace test {

Grid create_grid() {
    auto grid_config = eckit::YAMLConfiguration{ std::string{R"(
        type : "regional"
        nx : 35
        ny : 25
        north : -10
        south : -50
        east : 170
        west : 100
        y_numbering: +1
    )"}};
    return Grid{grid_config};
}
int truncation_x() {
    return (35-1)/2;
}
int truncation_y() {
    return (25-1)/2;
}

//-----------------------------------------------------------------------------

CASE("Test ectrans partitioner for LAM") {
    auto grid = create_grid();
    for( bool split_y : {false, true}) {
        SECTION("split_y = " + std::to_string(split_y)) {
            auto ectrans_partitioner       = grid::Partitioner("ectrans", util::Config("split_y",split_y));
            auto checkerboard_partitioner  = grid::Partitioner("checkerboard", util::Config("split_x",true)("split_y",split_y));
            auto ectrans_distribution      = grid::Distribution(grid, ectrans_partitioner);
            auto checkerboard_distribution = grid::Distribution(grid, checkerboard_partitioner);
            // Compared to checkerboard, ectrans has no control over split_x!
            EXPECT_EQ(ectrans_distribution.nb_pts(), checkerboard_distribution.nb_pts());
            EXPECT_EQ(ectrans_distribution.size(), checkerboard_distribution.size());
            int size = ectrans_distribution.size();
            for( int i=0; i<size; ++i) {
                EXPECT_EQ(ectrans_distribution.partition(i), checkerboard_distribution.partition(i));
            }

            functionspace::StructuredColumns fs_gp(grid, ectrans_partitioner);
            EXPECT_EQ(fs_gp.size(), ectrans_distribution.nb_pts()[mpi::rank()]);
        }
    }
}

#if 1
CASE("Test ectrans transform") {
    if (not trans::Trans::hasBackend("ectrans")) {
        Log::warning() << "Not testing as ectrans backend is not available" << std::endl;
        return;
    }
    auto grid = create_grid();
    bool split_y = false;
    trans::Trans trans(grid, truncation_x(), truncation_y(), option::split_y(split_y));

    functionspace::Spectral fs_sp(truncation_x(),truncation_y());
    functionspace::StructuredColumns fs_gp(grid, grid::Partitioner("ectrans", util::Config("split_y",split_y)));

    int computed_nspec2g = fs_sp.compute_nspec2g(truncation_x(),truncation_y());
    ATLAS_DEBUG_VAR(computed_nspec2g);
    int nspec2g = fs_sp.nb_spectral_coefficients_global();
    ATLAS_DEBUG_VAR(nspec2g);

    {
    int nvordiv = 1;
    int nscalar = 2;

    auto gp     = fs_gp.createField<double>(option::name("gp_scalar")|option::variables(nscalar));
    auto gpwind = fs_gp.createField<double>(option::name("gp_wind")  |option::variables(nvordiv*2));

    int mpi_root = 0;
    auto gpg = fs_gp.createField<double>(option::name("gp_scalar_global")|option::variables(nscalar)|option::global(mpi_root));
    auto gpwindg = fs_gp.createField<double>(option::name("gp_wind_global")|option::variables(nvordiv*2)|option::global(mpi_root));
    if (mpi::rank() == mpi_root) {
        auto rgpg = array::make_view<double,2>(gpg);
        auto rgpwindg = array::make_view<double,2>(gpwindg);
        ATLAS_DEBUG_VAR(rgpwindg.shape(0));
        ATLAS_DEBUG_VAR(rgpwindg.shape(1));
        for (int i=0; i<rgpwindg.shape(0); ++i) {
            for (int j=0; j<nvordiv; j++) {
                rgpwindg(i,2*j)   = 2*j;   // U
                rgpwindg(i,2*j+1) = 2*j+1; // V
            }
            for (int j=0; j<nscalar; j++) {
                rgpg(i,j) = 2*nvordiv+j; // scalar
            }
        }
    }
    fs_gp.scatter(gpg, gp);
    fs_gp.scatter(gpwindg, gpwind);

    if (mpi::rank() == mpi_root) {
        auto rgpg = array::make_view<double,2>(gpg);
        rgpg.assign(-1.);
    }
    
    auto sp    = fs_sp.createField<double>(option::name("sp_scalar")    |option::variables(nscalar));
    auto spvor = fs_sp.createField<double>(option::name("sp_vorticity") |option::variables(nvordiv));
    auto spdiv = fs_sp.createField<double>(option::name("sp_divergence")|option::variables(nvordiv));

    trans.dirtrans(gp, sp);
    trans.dirtrans_wind2vordiv(gpwind, spvor, spdiv);


    ATLAS_DEBUG_VAR(spvor.metadata().getDoubleVector("rmeanu"));
    ATLAS_DEBUG_VAR(spvor.metadata().getDoubleVector("rmeanv"));

    trans.invtrans(sp, gp);
    trans.invtrans_vordiv2wind(spvor, spdiv, gpwind);

    fs_gp.gather(gp, gpg);
    fs_gp.gather(gpwind, gpwindg);

    if (mpi::rank() == mpi_root) {
        auto rgpg = array::make_view<double,2>(gpg);
        rgpg.dump(std::cout);
    
        auto rgpwindg = array::make_view<double,2>(gpwindg);
        rgpwindg.dump(std::cout);
    }
    }

#if 0
{
  // Allocate gridpoint data
  int nvordiv = 1;
  int nscalar = 2;
  int nfld  = 2*nvordiv+nscalar;
  double* rgp  = malloc( sizeof(double) * nfld *trans.ngptot  );

  // Load data on proc 1
  double* rgpg = NULL;
  if( trans.myproc == 1 )
  {
    rgpg = malloc( sizeof(double) * 4*trans.ngptotg );
    int i, j;
    for ( j=0;j<nvordiv;j++) {
		  for( i=0; i<trans.ngptotg; ++i )
		  {
		    rgpg[(2*j)*trans.ngptotg+i] = 2*j; 		// U
		    rgpg[(2*j+1)*trans.ngptotg+i] = 2*j+1; // V
		  }
		}
    for ( j=0;j<nscalar;j++) {
		  for( i=0; i<trans.ngptotg; ++i )
		  {
		    rgpg[(2*nvordiv+j)*trans.ngptotg+i] =2*nvordiv+j; // scalar
		  }
		}
  }

  // Distribute gridpoint data from proc 1
  int* nfrom = malloc( sizeof(int) * nfld );
  for ( int j=0;j<nfld;j++) nfrom[j]=1;

  struct DistGrid_t distgrid = new_distgrid(&trans);
    distgrid.nfrom = nfrom;
    distgrid.rgpg  = rgpg;
    distgrid.rgp   = rgp;
    distgrid.nfld  = nfld;
  trans_distgrid(&distgrid);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
      if(nout < trans.ngptot-1) {
        printf("rgp[%d][...] : ...\n",j);
        i = trans.ngptot-1;
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
    }
  }

  // Allocate spectral data

  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
	double* rmeanu = malloc( sizeof(double) * nvordiv);
	double* rmeanv = malloc( sizeof(double) * nvordiv);
	
  // Direct Transform
  struct DirTrans_t dirtrans = new_dirtrans(&trans);
    dirtrans.nscalar   = nscalar;
    dirtrans.nvordiv   = nvordiv;
    dirtrans.rgp       = rgp;
    dirtrans.rspscalar = rspscalar;
    dirtrans.rspvor    = rspvor;
    dirtrans.rspdiv    = rspdiv;
    dirtrans.rmeanu    = rmeanu;
    dirtrans.rmeanv    = rmeanv;
  trans_dirtrans(&dirtrans);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalar[%d][%d] : %f\n",j,i,rspscalar[i*nscalar+j]);
    }
  }

  // Gather spectral field (for fun)
  int* nto = malloc( sizeof(int) * nscalar );
  nto[0] = 1;
  nto[1] = 1;

  double* rspscalarg = NULL;
  if( trans.myproc == 1 )
    rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );

  struct GathSpec_t gathspec = new_gathspec(&trans);
    gathspec.rspec  = rspscalar;
    gathspec.rspecg = rspscalarg;
    gathspec.nfld   = nscalar;
    gathspec.nto    = nto;
  trans_gathspec(&gathspec);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nvordiv; ++j)
    {
      printf("rmeanu[%d] : %f\n",j,rmeanu[j]);
      printf("rmeanv[%d] : %f\n",j,rmeanv[j]);
      for( i=0; i<nout; ++i )
        printf("rspvor[%d][%d] : %f\n",j,i,rspvor[i*nvordiv+j]);
        printf("rspdiv[%d][%d] : %f\n",j,i,rspdiv[i*nvordiv+j]);
    }
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalarg[%d][%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);
    }
  }
  
  // Distribute spectral field (for fun)
  struct DistSpec_t distspec = new_distspec(&trans);
    distspec.rspec  = rspscalar;
    distspec.rspecg = rspscalarg;
    distspec.nfld   = nscalar;
    distspec.nfrom  = nto;
  trans_distspec(&distspec);

  // Inverse Transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;
    invtrans.nvordiv   = nvordiv;
    invtrans.rspscalar = rspscalar;
    invtrans.rspvor    = rspvor;
    invtrans.rspdiv    = rspdiv;
    invtrans.rmeanu    = rmeanu;
    invtrans.rmeanv    = rmeanv;
    invtrans.rgp       = rgp;
  trans_invtrans(&invtrans);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
      if(nout < trans.ngptot-1) {
        printf("rgp[%d][...] : ...\n",j);
        i = trans.ngptot-1;
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
    }
  }

  // Gather gridpoint fields
  struct GathGrid_t gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nto  = nfrom;
    gathgrid.nfld = nfld;
  trans_gathgrid(&gathgrid);


  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgpg[%d][%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
      }
      if(nout < trans.ngptot-1) {
        printf("rgpg[%d][...] : ...\n",j);
        i = trans.ngptotg-1;
        printf("rgpg[%d][%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
      }
    }
  }

  // Deallocate arrays
  free(rgp);
  free(rgpg);
  free(rspscalar);
  free(rspscalarg);
  free(rspvor);
  free(rspdiv);
  free(nfrom);
  free(nto);
  free(rmeanu);
  free(rmeanv);
}

#endif



}
#endif

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
