/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "eckit/value/Properties.h"
#include "eckit/value/Params.h"
#include "eckit/memory/Owned.h"
#include "atlas/array/ArrayView.h"
#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
namespace field {
    class Field;
    class FieldSet;
}
}

namespace atlas {
namespace array {
    class Array;
}
}

namespace atlas {
namespace functionspace {
    class NodeColumns;
    class Spectral;
}
}

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

enum FFT { FFT992=TRANS_FFT992, FFTW=TRANS_FFTW };

//-----------------------------------------------------------------------------

class TransParameters : public eckit::Properties {
public:
  TransParameters() {}
  TransParameters(eckit::Stream& s) : eckit::Properties(s) { }
  ~TransParameters() {}
  static const char* className() { return "atlas::trans::TransParameters"; }
};

//-----------------------------------------------------------------------------

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
  Trans(const grid::Grid& g, const Options& = Options() );

  /// @brief Constructor given Gaussian N number for grid-only setup
  ///        This is equivalent to a (regular) Gaussian grid with N number
  ///        (e.g. for parallelisation routines)
  Trans( const size_t N, const Options& = Options() );

  /// @brief Constructor given grid and spectral truncation
  Trans( const grid::Grid& g, const size_t nsmax, const Options& = Options() );

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

  array::ArrayView<int,1> nloen() const
  {
    ASSERT( trans_.nloen != NULL );
    return array::ArrayView<int,1>(trans_.nloen, trans_.ndgl);
  }

  const int* n_regions(int& size) const
  {
    size = trans_.n_regions_NS;
    ASSERT( trans_.n_regions != NULL );
    return trans_.n_regions;
  }

  array::ArrayView<int,1> n_regions() const
  {
    ASSERT( trans_.n_regions != NULL );
    return array::ArrayView<int,1>(trans_.n_regions, trans_.n_regions_NS);
  }


  const int* nfrstlat(int& size) const
  {
    size = trans_.n_regions_NS;
    if( trans_.nfrstlat == NULL ) ::trans_inquire(&trans_,"nfrstlat");
    return trans_.nfrstlat;
  }

  array::ArrayView<int,1> nfrstlat() const
  {
    if( trans_.nfrstlat == NULL ) ::trans_inquire(&trans_,"nfrstlat");
    return array::ArrayView<int,1>(trans_.nfrstlat, trans_.n_regions_NS);
  }

  const int* nlstlat(int& size) const
  {
    size = trans_.n_regions_NS;
    if( trans_.nlstlat == NULL ) ::trans_inquire(&trans_,"nlstlat");
    return trans_.nlstlat;
  }

  array::ArrayView<int,1> nlstlat() const
  {
    if( trans_.nlstlat == NULL ) ::trans_inquire(&trans_,"nlstlat");
    return array::ArrayView<int,1>(trans_.nlstlat, trans_.n_regions_NS);
  }

  const int* nptrfrstlat(int& size) const
  {
    size = trans_.n_regions_NS;
    if( trans_.nptrfrstlat == NULL ) ::trans_inquire(&trans_,"nptrfrstlat");
    return trans_.nptrfrstlat;
  }

  array::ArrayView<int,1> nptrfrstlat() const
  {
    if( trans_.nptrfrstlat == NULL ) ::trans_inquire(&trans_,"nptrfrstlat");
    return array::ArrayView<int,1>(trans_.nptrfrstlat, trans_.n_regions_NS);
  }

  const int* nsta(int& sizef2, int& sizef1) const
  {
    sizef1 = trans_.ndgl+trans_.n_regions_NS-1;
    sizef2 = trans_.n_regions_EW;
    if( trans_.nsta == NULL ) ::trans_inquire(&trans_,"nsta");
    return trans_.nsta;
  }

  array::ArrayView<int,2> nsta() const
  {
    if( trans_.nsta == NULL ) ::trans_inquire(&trans_,"nsta");
    return array::ArrayView<int,2>( trans_.nsta, array::make_shape(trans_.n_regions_EW, trans_.ndgl+trans_.n_regions_NS-1) );
  }

  const int* nonl(int& sizef2, int& sizef1) const
  {
    sizef1 = trans_.ndgl+trans_.n_regions_NS-1;
    sizef2 = trans_.n_regions_EW;
    if( trans_.nonl == NULL ) ::trans_inquire(&trans_,"nonl");
    return trans_.nonl;
  }

  array::ArrayView<int,2> nonl() const
  {
    if( trans_.nonl == NULL ) ::trans_inquire(&trans_,"nonl");
    return array::ArrayView<int,2>( trans_.nonl, array::make_shape(trans_.n_regions_EW, trans_.ndgl+trans_.n_regions_NS-1) );
  }

  const int* nmyms(int& size) const
  {
    size = trans_.nump;
    if( trans_.nmyms == NULL ) ::trans_inquire(&trans_,"nmyms");
    return trans_.nmyms;
  }

  array::ArrayView<int,1> nmyms() const
  {
    if( trans_.nmyms == NULL ) ::trans_inquire(&trans_,"nmyms");
    return array::ArrayView<int,1> (trans_.nmyms, trans_.nump);
  }

  const int* nasm0(int& size) const
  {
    size = trans_.nsmax+1; // +1 because zeroth wave included
    if( trans_.nasm0 == NULL ) ::trans_inquire(&trans_,"nasm0");
    return trans_.nasm0;
  }

  array::ArrayView<int,1> nasm0() const
  {
    if( trans_.nasm0 == NULL ) ::trans_inquire(&trans_,"nasm0");
    return array::ArrayView<int,1> (trans_.nasm0, trans_.nsmax+1);
  }

  const int* nvalue(int& size) const
  {
    size = trans_.nspec2;
    if( trans_.nvalue == NULL ) ::trans_inquire(&trans_,"nvalue");
    return trans_.nvalue;
  }

  array::ArrayView<int,1> nvalue() const
  {
    if( trans_.nvalue == NULL ) ::trans_inquire(&trans_,"nvalue");
    return array::ArrayView<int,1> (trans_.nvalue, trans_.nspec2);
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

  void dirtrans(const field::Field& gpfield,
                      field::Field& spfield,
                const TransParameters& = TransParameters()) const;
  void dirtrans(const field::FieldSet& gpfields,
                      field::FieldSet& spfields,
                const TransParameters& = TransParameters()) const;

  void dirtrans(const functionspace::NodeColumns&,    const field::Field& gpfield,
                const functionspace::Spectral&,       field::Field& spfield,
                const TransParameters& = TransParameters()) const;
  void dirtrans(const functionspace::NodeColumns&,    const field::FieldSet& gpfields,
                const functionspace::Spectral&,       field::FieldSet& spfields,
                const TransParameters& = TransParameters()) const;
  void dirtrans_wind2vordiv(const functionspace::NodeColumns&, const field::Field& gpwind,
                            const functionspace::Spectral&, field::Field& spvor, field::Field& spdiv,
                            const TransParameters& = TransParameters()) const;

  void invtrans(const field::Field& spfield,
                      field::Field& gpfield,
                const TransParameters& = TransParameters()) const;
  void invtrans(const field::FieldSet& spfields,
                      field::FieldSet& gpfields,
                const TransParameters& = TransParameters()) const;

  void invtrans(const functionspace::Spectral&, const field::Field& spfield,
                const functionspace::NodeColumns&,          field::Field& gpfield,
                const TransParameters& = TransParameters()) const;
  void invtrans(const functionspace::Spectral&, const field::FieldSet& spfields,
                const functionspace::NodeColumns&,          field::FieldSet& gpfields,
                const TransParameters& = TransParameters()) const;
  void invtrans_vordiv2wind(const functionspace::Spectral&, const field::Field& spvor, const field::Field& spdiv,
                            const functionspace::NodeColumns&, field::Field& gpwind,
                            const TransParameters& = TransParameters()) const;


private:

  void ctor_rgg(const size_t nlat, const long pl[], size_t nsmax, const Options& );

  void ctor_lonlat(const size_t nlon, const size_t nlat, size_t nsmax, const Options& );


private:
  mutable Trans_t trans_;
};

//-----------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
#define functionspace_NodeColumns functionspace::NodeColumns
#define functionspace_Spectral functionspace::Spectral
#define grid_Grid grid::Grid::grid_t
#define field_Field field::Field
#define field_FieldSet field::FieldSet

extern "C"
{
  Trans* atlas__Trans__new (const grid_Grid* grid, int nsmax);
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
  void atlas__Trans__dirtrans_fieldset (const Trans* This, const field_FieldSet* gpfields, field_FieldSet* spfields, const TransParameters* parameters);
  void atlas__Trans__dirtrans_field (const Trans* This, const field_Field* gpfield, field_Field* spfield, const TransParameters* parameters);
  void atlas__Trans__invtrans_fieldset (const Trans* This, const field_FieldSet* spfields, field_FieldSet* gpfields, const TransParameters* parameters);
  void atlas__Trans__invtrans_field (const Trans* This, const field_Field* spfield, field_Field* gpfield, const TransParameters* parameters);
  void atlas__Trans__dirtrans_fieldset_nodes (const Trans* This, const functionspace_NodeColumns* gp, const field_FieldSet* gpfields, const functionspace_Spectral* sp, field_FieldSet* spfields, const TransParameters* parameters);
  void atlas__Trans__invtrans_fieldset_nodes (const Trans* This, const functionspace_Spectral* sp, const field_FieldSet* spfields, const functionspace_NodeColumns* gp, field_FieldSet* gpfields, const TransParameters* parameters);
  void atlas__Trans__dirtrans_field_nodes (const Trans* This, const functionspace_NodeColumns* gp, const field_Field* gpfield, const functionspace_Spectral* sp, field_Field* spfield, const TransParameters* parameters);
  void atlas__Trans__invtrans_field_nodes (const Trans* This, const functionspace_Spectral* sp, const field_Field* spfield, const functionspace_NodeColumns* gp, field_Field* gpfield, const TransParameters* parameters);
  void atlas__Trans__dirtrans_wind2vordiv_field_nodes (const Trans* This, const functionspace_NodeColumns* gp, const field_Field* gpwind, const functionspace_Spectral* sp, field_Field* spvor, field_Field* spdiv, const TransParameters* parameters);
  void atlas__Trans__invtrans_vordiv2wind_field_nodes (const Trans* This, const functionspace_Spectral* sp, const field_Field* spvor, const field_Field* spdiv, const functionspace_NodeColumns* gp, field_Field* gpwind, const TransParameters* parameters);
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
#undef grid_Grid
#undef field_FieldSet
#undef field_Field
#undef functionspace_NodeColumns
#undef functionspace_Spectral

// ------------------------------------------------------------------

} // namespace trans
} // namespace atlas


#endif // atlas_trans_Trans_h
