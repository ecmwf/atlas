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

#include "eckit/value/Properties.h"
#include "eckit/value/Params.h"
#include "transi/trans.h"

namespace atlas {
namespace grids {
  class ReducedGrid;
}
}

namespace atlas {
namespace trans {

enum FFT { FFT992=TRANS_FFT992, FFTW=TRANS_FFTW };

class Trans {
private:
  typedef struct ::Trans_t Trans_t;

public:

  class Options : public eckit::Properties {
  public:

    Options();
    Options(eckit::Stream& s) : eckit::Properties(s) {}
    ~Options() {}
    static const char* className() { return "atlas::trans::Trans::Options"; }

    void set_split_latitudes(bool);
    void set_fft( FFT );
    void set_flt(bool);

    bool split_latitudes() const;
    FFT fft() const;
    bool flt() const;

    friend std::ostream& operator<<( std::ostream& os, const Options& p) { p.print(os); return os;}

  private:
    friend eckit::Params::value_t get( const Options& p, const eckit::Params::key_t& key );
    friend void print( const Options& p, std::ostream& s );
    friend void encode( const Options& p, eckit::Stream& s );

    void print( std::ostream& ) const;
  };

public:

  /// @brief Constructor for grid-only setup
  ///        (e.g. for parallelisation routines)
  Trans(const grids::ReducedGrid& g, const Options& = Options() );

  /// @brief Constructor given Gaussian N number for grid-only setup
  ///        This is equivalent to a (regular) Gaussian grid with N number
  ///        (e.g. for parallelisation routines)
  Trans( const int N, const Options& = Options() );

  /// @brief Constructor given grid and spectral truncation
  Trans( const grids::ReducedGrid& g, const int nsmax, const Options& = Options() );

  /// @brief Constructor given npts_per_lat (aka PL array) and spectral truncation
  Trans( const std::vector<int>& npts_per_lat, const int nsmax, const Options& = Options() );

  /// @brief Constructor given Gaussian N number and spectral truncation
  ///        This is equivalent to a (regular) Gaussian grid with N number
  Trans( const int N, const int nsmax, const Options& = Options() );

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

  void ctor_rgg(const int ndgl, const int nloen[], int nsmax, const Options& );

  void ctor_lonlat(const int nlon, const int nlat, int nsmax, const Options& );


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
