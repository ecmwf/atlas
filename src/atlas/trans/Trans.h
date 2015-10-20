/*
 * (C) Copyright 1996-2015 ECMWF.
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
#include "eckit/memory/Owned.h"
#include "transi/trans.h"
#include "atlas/util/ArrayView.h"

namespace atlas {
  class FieldSet;
  class Grid;
  class Field;
  class Array;
  namespace functionspace {
    class Nodes;
    class Spectral;
  }
}

namespace atlas {
namespace trans {

enum FFT { FFT992=TRANS_FFT992, FFTW=TRANS_FFTW };

class TransParameters : public eckit::Properties {
public:
  TransParameters() {}
  TransParameters(eckit::Stream& s) : eckit::Properties(s) { }
  ~TransParameters() {}
  static const char* className() { return "atlas::trans::TransParameters"; }
};


class Trans : public eckit::Owned {
private:
  typedef struct ::Trans_t Trans_t;

public:

  class Options : public eckit::Properties {
  public:

    Options();
    Options(eckit::Stream& s) : eckit::Properties(s) { set_cache(NULL,0); }
    ~Options() {}
    static const char* className() { return "atlas::trans::Trans::Options"; }

    void set_split_latitudes(bool);
    void set_fft( FFT );
    void set_flt(bool);
    void set_cache(const void* buffer, const size_t size);
    void set_read(const std::string&);
    void set_write(const std::string&);

    bool split_latitudes() const;
    FFT fft() const;
    bool flt() const;
    const void* cache()const;
    size_t cachesize() const;
    std::string read() const;
    std::string write() const;

    friend std::ostream& operator<<( std::ostream& os, const Options& p) { p.print(os); return os;}

  private:
    friend eckit::Params::value_t getValue( const Options& p, const eckit::Params::key_t& key );
    friend void print( const Options& p, std::ostream& s );
    friend void encode( const Options& p, eckit::Stream& s );

    void print( std::ostream& ) const;

  private: // not represented in Property internals
    const void*  cacheptr_;
    size_t cachesize_;
  };

public:

  /// @brief Constructor for grid-only setup
  ///        (e.g. for parallelisation routines)
  Trans(const Grid& g, const Options& = Options() );

  /// @brief Constructor given Gaussian N number for grid-only setup
  ///        This is equivalent to a (regular) Gaussian grid with N number
  ///        (e.g. for parallelisation routines)
  Trans( const size_t N, const Options& = Options() );

  /// @brief Constructor given grid and spectral truncation
  Trans( const Grid& g, const size_t nsmax, const Options& = Options() );

  /// @brief Constructor given npts_per_lat (aka PL array) and spectral truncation
  Trans( const std::vector<size_t>& npts_per_lat, const size_t nsmax, const Options& = Options() );

  /// @brief Constructor given Gaussian N number and spectral truncation
  ///        This is equivalent to a (regular) Gaussian grid with N number
  Trans( const size_t N, const size_t nsmax, const Options& = Options() );

  virtual ~Trans();
  operator Trans_t*() const { return &trans_; }

  int        handle()       const { return trans_.handle; }
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
  int        nump()         const { return trans_.nump; }
  int        nproc()        const { return trans_.nproc; }
  int        myproc(int proc0=0) const { return trans_.myproc-1+proc0; }

  const int* nloen(int& size) const
  {
    size = trans_.ndgl;
    ASSERT( trans_.nloen != NULL );
    return trans_.nloen;
  }

  ArrayView<int,1> nloen() const
  {
    ASSERT( trans_.nloen != NULL );
    return ArrayView<int,1>(trans_.nloen, trans_.ndgl);
  }

  const int* n_regions(int& size) const
  {
    size = trans_.n_regions_NS;
    ASSERT( trans_.n_regions != NULL );
    return trans_.n_regions;
  }

  ArrayView<int,1> n_regions() const
  {
    ASSERT( trans_.n_regions != NULL );
    return ArrayView<int,1>(trans_.n_regions, trans_.n_regions_NS);
  }


  const int* nfrstlat(int& size) const
  {
    size = trans_.n_regions_NS;
    if( trans_.nfrstlat == NULL ) ::trans_inquire(&trans_,"nfrstlat");
    return trans_.nfrstlat;
  }

  ArrayView<int,1> nfrstlat() const
  {
    if( trans_.nfrstlat == NULL ) ::trans_inquire(&trans_,"nfrstlat");
    return ArrayView<int,1>(trans_.nfrstlat, trans_.n_regions_NS);
  }

  const int* nlstlat(int& size) const
  {
    size = trans_.n_regions_NS;
    if( trans_.nlstlat == NULL ) ::trans_inquire(&trans_,"nlstlat");
    return trans_.nlstlat;
  }

  ArrayView<int,1> nlstlat() const
  {
    if( trans_.nlstlat == NULL ) ::trans_inquire(&trans_,"nlstlat");
    return ArrayView<int,1>(trans_.nlstlat, trans_.n_regions_NS);
  }

  const int* nptrfrstlat(int& size) const
  {
    size = trans_.n_regions_NS;
    if( trans_.nptrfrstlat == NULL ) ::trans_inquire(&trans_,"nptrfrstlat");
    return trans_.nptrfrstlat;
  }

  ArrayView<int,1> nptrfrstlat() const
  {
    if( trans_.nptrfrstlat == NULL ) ::trans_inquire(&trans_,"nptrfrstlat");
    return ArrayView<int,1>(trans_.nptrfrstlat, trans_.n_regions_NS);
  }

  const int* nsta(int& sizef2, int& sizef1) const
  {
    sizef1 = trans_.ndgl+trans_.n_regions_NS-1;
    sizef2 = trans_.n_regions_EW;
    if( trans_.nsta == NULL ) ::trans_inquire(&trans_,"nsta");
    return trans_.nsta;
  }

  ArrayView<int,2> nsta() const
  {
    if( trans_.nsta == NULL ) ::trans_inquire(&trans_,"nsta");
    return ArrayView<int,2>( trans_.nsta, make_shape(trans_.n_regions_EW, trans_.ndgl+trans_.n_regions_NS-1) );
  }

  const int* nonl(int& sizef2, int& sizef1) const
  {
    sizef1 = trans_.ndgl+trans_.n_regions_NS-1;
    sizef2 = trans_.n_regions_EW;
    if( trans_.nonl == NULL ) ::trans_inquire(&trans_,"nonl");
    return trans_.nonl;
  }

  ArrayView<int,2> nonl() const
  {
    if( trans_.nonl == NULL ) ::trans_inquire(&trans_,"nonl");
    return ArrayView<int,2>( trans_.nonl, make_shape(trans_.n_regions_EW, trans_.ndgl+trans_.n_regions_NS-1) );
  }

  const int* nmyms(int& size) const
  {
    size = trans_.nump;
    if( trans_.nmyms == NULL ) ::trans_inquire(&trans_,"nmyms");
    return trans_.nmyms;
  }

  ArrayView<int,1> nmyms() const
  {
    if( trans_.nmyms == NULL ) ::trans_inquire(&trans_,"nmyms");
    return ArrayView<int,1> (trans_.nmyms, trans_.nump);
  }

  const int* nasm0(int& size) const
  {
    size = trans_.nsmax+1; // +1 because zeroth wave included
    if( trans_.nasm0 == NULL ) ::trans_inquire(&trans_,"nasm0");
    return trans_.nasm0;
  }

  ArrayView<int,1> nasm0() const
  {
    if( trans_.nasm0 == NULL ) ::trans_inquire(&trans_,"nasm0");
    return ArrayView<int,1> (trans_.nasm0, trans_.nsmax+1);
  }

  const int* nvalue(int& size) const
  {
    size = trans_.nspec2;
    if( trans_.nvalue == NULL ) ::trans_inquire(&trans_,"nvalue");
    return trans_.nvalue;
  }

  ArrayView<int,1> nvalue() const
  {
    if( trans_.nvalue == NULL ) ::trans_inquire(&trans_,"nvalue");
    return ArrayView<int,1> (trans_.nvalue, trans_.nspec2);
  }

public:


  /*!
   * @brief distspec
   * @param nb_fields
   * @param origin
   * @param global_spectra
   * @param spectra
   */
  void distspec( const int nb_fields, const int origin[], const double global_spectra[], double spectra[] ) const;

  /*!
   * @brief gathspec
   * @param nb_fields
   * @param destination
   * @param spectra
   * @param global_spectra
   */
  void gathspec( const int nb_fields, const int destination[], const double spectra[], double global_spectra[] ) const;

  /*!
   * @brief distgrid
   * @param nb_fields
   * @param origin
   * @param global_fields
   * @param fields
   */
  void distgrid( const int nb_fields, const int origin[], const double global_fields[], double fields[] ) const;

  /*!
   * @brief gathgrid
   * @param nb_fields
   * @param destination
   * @param fields
   * @param global_fields
   */
  void gathgrid( const int nb_fields, const int destination[], const double fields[], double global_fields[] ) const;

  /*!
   * @brief invtrans
   * @param nb_fields
   * @param scalar_spectra
   * @param scalar_fields
   */
  void invtrans( const int nb_fields, const double scalar_spectra[], double scalar_fields[] ) const;

  /*!
   * @brief Inverse transform of vorticity/divergence to wind(U/V)
   * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
   */
  void invtrans( const int nb_fields, const double vorticity_spectra[], const double divergence_spectra[], double wind_fields[] ) const;

  /*!
   * @brief Direct transform of scalar fields
   */
  void dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[] ) const;

  /*!
   * @brief Direct transform of wind(U/V) to vorticity/divergence
   * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
   */
  void dirtrans(const int nb_fields, const double wind_fields[], double vorticity_spectra[], double divergence_spectra[] ) const;

  void dirtrans(const functionspace::Nodes&,    const Field& gpfield,
                const functionspace::Spectral&,       Field& spfield,
                const TransParameters& = TransParameters()) const;
  void dirtrans(const functionspace::Nodes&,    const FieldSet& gpfields,
                const functionspace::Spectral&,       FieldSet& spfields,
                const TransParameters& = TransParameters()) const;
  void dirtrans_wind2vordiv(const functionspace::Nodes&, const Field& gpwind,
                            const functionspace::Spectral&, Field& spvor, Field& spdiv,
                            const TransParameters& = TransParameters()) const;

  void invtrans(const functionspace::Spectral&, const Field& spfield,
                const functionspace::Nodes&,          Field& gpfield,
                const TransParameters& = TransParameters()) const;
  void invtrans(const functionspace::Spectral&, const FieldSet& spfields,
                const functionspace::Nodes&,          FieldSet& gpfields,
                const TransParameters& = TransParameters()) const;
  void invtrans_vordiv2wind(const functionspace::Spectral&, const Field& spvor, const Field& spdiv,
                            const functionspace::Nodes&, Field& gpwind,
                            const TransParameters& = TransParameters()) const;


private:

  void ctor_rgg(const size_t ndgl, const int nloen[], size_t nsmax, const Options& );

  void ctor_lonlat(const size_t nlon, const size_t nlat, size_t nsmax, const Options& );


private:
  mutable Trans_t trans_;
};



// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

#define NodesFunctionSpace functionspace::Nodes
#define SpectralFunctionSpace functionspace::Spectral
extern "C"
{
  Trans* atlas__Trans__new (const Grid* grid, int nsmax);
  void atlas__Trans__delete (Trans* trans);
  int atlas__Trans__handle (const Trans* trans);
  void atlas__Trans__distspec (const Trans* t, int nb_fields, int origin[], double global_spectra[], double spectra[]);
  void atlas__Trans__gathspec (const Trans* t, int nb_fields, int destination[], double spectra[], double global_spectra[]);
  void atlas__Trans__distgrid (const Trans* t, int nb_fields, int origin[], double global_fields[], double fields[]);
  void atlas__Trans__gathgrid (const Trans* t, int nb_fields, int destination[], double fields[], double global_fields[]);
  void atlas__Trans__invtrans_scalar (const Trans* t, int nb_fields, double scalar_spectra[], double scalar_fields[]);
  void atlas__Trans__invtrans_vordiv2wind (const Trans* t, int nb_fields, double vorticity_spectra[], double divergence_spectra[], double wind_fields[]);
  void atlas__Trans__dirtrans_scalar (const Trans* t, int nb_fields, double scalar_fields[], double scalar_spectra[]);
  void atlas__Trans__dirtrans_wind2vordiv (const Trans* t, int nb_fields, double wind_fields[], double vorticity_spectra[], double divergence_spectra[]);
  void atlas__Trans__dirtrans_fieldset (const Trans* This, const NodesFunctionSpace* gp, const FieldSet* gpfields, const SpectralFunctionSpace* sp, FieldSet* spfields, const TransParameters* parameters);
  void atlas__Trans__invtrans_fieldset (const Trans* This, const SpectralFunctionSpace* sp, const FieldSet* spfields, const NodesFunctionSpace* gp, FieldSet* gpfields, const TransParameters* parameters);
  void atlas__Trans__dirtrans_field (const Trans* This, const NodesFunctionSpace* gp, const Field* gpfield, const SpectralFunctionSpace* sp, Field* spfield, const TransParameters* parameters);
  void atlas__Trans__invtrans_field (const Trans* This, const SpectralFunctionSpace* sp, const Field* spfield, const NodesFunctionSpace* gp, Field* gpfield, const TransParameters* parameters);
  void atlas__Trans__dirtrans_wind2vordiv_field (const Trans* This, const NodesFunctionSpace* gp, const Field* gpwind, const SpectralFunctionSpace* sp, Field* spvor, Field* spdiv, const TransParameters* parameters);
  void atlas__Trans__invtrans_vordiv2wind_field (const Trans* This, const SpectralFunctionSpace* sp, const Field* spvor, const Field* spdiv, const NodesFunctionSpace* gp, Field* gpwind, const TransParameters* parameters);
  int atlas__Trans__nproc (const Trans* This);
  int atlas__Trans__myproc (const Trans* This, int proc0);
  int atlas__Trans__ndgl (const Trans* This);
  int atlas__Trans__nsmax (const Trans* This);
  int atlas__Trans__ngptot (const Trans* This);
  int atlas__Trans__ngptotg (const Trans* This);
  int atlas__Trans__ngptotmx (const Trans* This);
  int atlas__Trans__nspec (const Trans* This);
  int atlas__Trans__nspec2 (const Trans* This);
  int atlas__Trans__nspec2g (const Trans* This);
  int atlas__Trans__nspec2mx (const Trans* This);
  int atlas__Trans__n_regions_NS (const Trans* This);
  int atlas__Trans__n_regions_EW (const Trans* This);
  int atlas__Trans__nump (const Trans* This);
  const int* atlas__Trans__nloen (const Trans* This, int &size);
  const int* atlas__Trans__n_regions (const Trans* This, int &size);
  const int* atlas__Trans__nfrstlat (const Trans* This, int &size);
  const int* atlas__Trans__nlstlat (const Trans* This, int &size);
  const int* atlas__Trans__nptrfrstlat (const Trans* This, int &size);
  const int* atlas__Trans__nsta (const Trans* This, int &sizef2, int &sizef1);
  const int* atlas__Trans__nonl (const Trans* This, int &sizef2, int &sizef1);
  const int* atlas__Trans__nmyms (const Trans* This, int &size);
  const int* atlas__Trans__nasm0 (const Trans* This, int &size);
  const int* atlas__Trans__nvalue (const Trans* This, int &size);
  TransParameters* atlas__TransParameters__new ();
  void atlas__TransParameters__delete (TransParameters* parameters);
}
#undef NodesFunctionSpace
#undef SpectralFunctionSpace

// ------------------------------------------------------------------


}
}


#endif // atlas_trans_Trans_h
