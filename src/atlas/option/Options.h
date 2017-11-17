/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/util/Config.h"
#include "atlas/util/Earth.h"
#include "atlas/array/DataType.h"

// ----------------------------------------------------------------------------

namespace atlas {
namespace option {

// ----------------------------------------------------------------------------

class type : public util::Config {
public:
  type( const std::string& );
};

// ----------------------------------------------------------------------------

class global : public util::Config {
public:
  global( size_t owner = 0 );
};

// ----------------------------------------------------------------------------

class levels : public util::Config {
public:
  levels( size_t );
};

// ----------------------------------------------------------------------------

class variables : public util::Config {
public:
  variables( size_t );
};

// ----------------------------------------------------------------------------

class name : public util::Config {
public:
  name( const std::string& );
};

// ----------------------------------------------------------------------------

template< typename T >
class datatypeT : public util::Config {
public:
  datatypeT();
};

// ----------------------------------------------------------------------------

class datatype : public util::Config {
public:
  datatype( array::DataType::kind_t );
  datatype( const std::string& );
  datatype( array::DataType );
};

// ----------------------------------------------------------------------------

class halo : public util::Config {
public:
  halo(size_t size);
};

// ----------------------------------------------------------------------------

class radius : public util::Config {
public:
  radius( double );
  radius( const std::string& = "Earth" );
};

// ----------------------------------------------------------------------------

class scalar_derivatives : public util::Config {
public:
  scalar_derivatives( bool );
};

// ----------------------------------------------------------------------------

class wind_EW_derivatives : public util::Config {
public:
  wind_EW_derivatives( bool );
};

// ----------------------------------------------------------------------------

class vorticity_divergence_fields : public util::Config {
public:
  vorticity_divergence_fields( bool );
};

// ----------------------------------------------------------------------------

class flt : public util::Config {
public:
  flt( bool );
};

// ----------------------------------------------------------------------------

enum FFT { FFT992=1, FFTW=2 };

class fft : public util::Config {
public:
  fft( FFT );
  fft( const std::string& );
};

// ----------------------------------------------------------------------------

class split_latitudes : public util::Config {
public:
  split_latitudes( bool );
};

// ----------------------------------------------------------------------------

class write_legendre : public util::Config {
public:
  write_legendre( const eckit::PathName& );
};

// ----------------------------------------------------------------------------

class read_legendre : public util::Config {
public:
  read_legendre( const eckit::PathName& );
};

// ----------------------------------------------------------------------------
// Definitions
// ----------------------------------------------------------------------------

template<typename T>
datatypeT<T>::datatypeT() {
  set("datatype",array::DataType::kind<T>());
}



} // namespace option
} // namespace atlas
