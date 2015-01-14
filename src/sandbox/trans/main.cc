#include <limits>
#include <cassert>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>
#include <stdlib.h> /* malloc, free */

#include "trans_api.h"

//------------------------------------------------------------------------------------------------------

/*! @example trans_sptogp.c 
 *
 * Transform spectral to gridpoint
 * 
 * This is an example of how to setup and
 * transform global spectral data to 
 * global gridpoint data
 */

// Following dummy functions are implementation details 
// that don't contribute to this example. They could be
// replaced with grib_api functionality
void read_grid( Trans* trans );
void read_rspecg( Trans* trans, double* rspecg[], int* nfrom[], int* nfld );
void write_rgpg( Trans* trans, double* rgpg[], int nfld );


int main ( int arc, char **argv ) 
{
  int jfld;
  Trans trans;
 
  // Read resolution information
  read_grid(&trans);

  // Register resolution in trans library
  trans_setup(&trans);

  // Declare global spectral data
  int nfld;
  double* rspecg = NULL;
  int* nfrom = NULL;

  // Read global spectral data (could be from grib file)
  read_rspecg(&trans,&rspecg,&nfrom,&nfld);

  // Distribute data to all procs
  double* rspec  = (double*) malloc( sizeof(double) * nfld *trans.nspec2  );
  DistSpec distspec = new_distspec(&trans);
    distspec.nfrom  = nfrom;
    distspec.rspecg = rspecg;
    distspec.rspec  = rspec;
    distspec.nfld   = nfld;
  trans_distspec(&distspec);


  // Transform sp to gp fields
  double* rgp = (double*) malloc( sizeof(double) * nfld*trans.ngptot );

  InvTrans invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nfld;
    invtrans.rspscalar = rspec;
    invtrans.rgp       = rgp;
  trans_invtrans(&invtrans);

  
  // Gather all gridpoint fields
  double* rgpg = NULL;
  if( trans.myproc == 1 )
    rgpg = (double*) malloc( sizeof(double) * nfld*trans.ngptotg );

  int* nto =  (int*) malloc( sizeof(int) * nfld );
  for( jfld=0; jfld<nfld; ++jfld )
    nto[jfld] = 1;

  GathGrid gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nfld = nfld;
    gathgrid.nto  = nto;
  trans_gathgrid(&gathgrid);

  // Write global spectral data (could be to grib file)
  if( trans.myproc == 1 )
    write_rgpg(&trans,&rgpg,nfld);

  // Deallocate and finalize
  free(rgp);
  free(rgpg);
  free(rspec);
  free(rspecg);
  free(nfrom);
  free(nto);

  trans_delete(&trans);

  trans_finalize();

  return 0;
}

//---------------------------------------------------------------------------
// Dummy functions, used in this example

void read_grid(Trans* trans)
{
  int i;
  int T159[] = {
     18,  25,  36,  40,  45,  50,  60,  64,  72,  72,
     80,  90,  96, 100, 108, 120, 120, 125, 135, 144,
    144, 150, 160, 160, 180, 180, 180, 192, 192, 200,
    200, 216, 216, 216, 225, 225, 240, 240, 240, 243,
    250, 250, 256, 270, 270, 270, 288, 288, 288, 288,
    288, 288, 300, 300, 300, 300, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 300, 300, 300, 300, 288, 288,
    288, 288, 288, 288, 270, 270, 270, 256, 250, 250,
    243, 240, 240, 240, 225, 225, 216, 216, 216, 200,
    200, 192, 192, 180, 180, 180, 160, 160, 150, 144,
    144, 135, 125, 120, 120, 108, 100,  96,  90,  80,
     72,  72,  64,  60,  50,  45,  40,  36,  25,  18,
  };
  trans->ndgl  = sizeof(T159)/sizeof(int);
  trans->nloen = (int*) malloc( sizeof(T159) );
  for( i=0; i<trans->ndgl; i++)  trans->nloen[i] = T159[i];

  // Assume Linear Grid
  trans->nsmax=(2*trans->ndgl-1)/2;
}


void read_rspecg(Trans* trans, double* rspecg[], int* nfrom[], int* nfld )
{
  int i;
  int jfld;
  if( trans->myproc == 1 ) printf("read_rspecg ...\n");
  *nfld = 2;
  if( trans->myproc == 1 )
  {
    *rspecg = (double*) malloc( sizeof(double) * (*nfld)*trans->nspec2g );
    for( i=0; i<trans->nspec2g; ++i )
    {
      (*rspecg)[i*(*nfld) + 0] = (i==0 ? 1. : 0.); // scalar field 1
      (*rspecg)[i*(*nfld) + 1] = (i==0 ? 2. : 0.); // scalar field 2
    } 
  }
  *nfrom = (int*) malloc( sizeof(int) * (*nfld) );
  for (jfld=0; jfld<(*nfld); ++jfld)
    (*nfrom)[jfld] = 1;
  if( trans->myproc == 1 ) printf("read_rspecg ... done\n");
}

void write_rgpg(Trans* trans, double* rgpg[], int nfld )
{
  int jfld;
  if( trans->myproc == 1 ) printf("write_rgpg ...\n");
  for( jfld=0; jfld<nfld; ++jfld )
  {
    // output global field rgpg[jfld]
  }
  if( trans->myproc == 1 ) printf("write_rgpg ... done\n");
}

