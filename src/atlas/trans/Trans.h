/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_trans_Trans_h
#define atlas_trans_Trans_h

#include "transi/trans.h"

namespace atlas {
namespace grids {
  class ReducedGrid;
}
}

namespace atlas {
namespace trans {

class Trans {
private:
  typedef struct ::Trans_t Trans_t;

public:
  Trans(const grids::ReducedGrid& g);
  Trans(const grids::ReducedGrid& g, const int nsmax );
  Trans( const std::vector<int>& npts_per_lat, const int nsmax );
  virtual ~Trans();
  operator Trans_t*() const { return &trans_; }

  int        handle()       const { return trans_.handle; }
  int        nproc()        const { return trans_.nproc; }
  int        myproc()       const { return trans_.myproc; }
  int        ndgl()         const { return trans_.ndgl; }
  int        nsmax()        const { return trans_.nsmax; }
  int        ngptot()       const { return trans_.ngptot; }
  int        ngptotg()      const { return trans_.ngptotg; }
  int        ngptotmx()     const { return trans_.ngptotmx; }
  int        nspec()        const { return trans_.nspec; }
  int        nspec2()       const { return trans_.nspec2; }
  int        nspec2g()      const { return trans_.nspec2g; }
  int        nspec2mx()     const { return trans_.nspec2mx; }
  int        n_regions_NS() const { return trans_.n_regions_NS; }
  int        n_regions_EW() const { return trans_.n_regions_EW; }
  const int* nloen()        const { ASSERT( trans_.nloen     != NULL ); return trans_.nloen; }
  const int* n_regions()    const { ASSERT( trans_.n_regions != NULL ); return trans_.n_regions; }
  const int* nfrstlat()     const { if( trans_.nfrstlat    == NULL ) ::trans_inquire(&trans_,"nfrstlat");    return trans_.nfrstlat; }
  const int* nlstlat()      const { if( trans_.nlstlat     == NULL ) ::trans_inquire(&trans_,"nlstlat");     return trans_.nlstlat; }
  const int* nptrfrstlat()  const { if( trans_.nptrfrstlat == NULL ) ::trans_inquire(&trans_,"nptrfrstlat"); return trans_.nptrfrstlat; }
  const int* nsta()         const { if( trans_.nsta        == NULL ) ::trans_inquire(&trans_,"nsta");        return trans_.nsta; }
  const int* nonl()         const { if( trans_.nonl        == NULL ) ::trans_inquire(&trans_,"nonl");        return trans_.nonl; }

private:

  void ctor(const int ndgl, const int nloen[], int nsmax);

private:
  mutable Trans_t trans_;
};

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

typedef grids::ReducedGrid grids__ReducedGrid;

extern "C"
{
  Trans* atlas__Trans__new (grids__ReducedGrid* grid);
  void atlas__Trans__delete (Trans* trans);
  int atlas__Trans__handle (Trans* trans);
}
// ------------------------------------------------------------------


}
}


#endif // atlas_trans_Trans_h
