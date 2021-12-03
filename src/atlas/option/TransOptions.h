/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/util/Config.h"

namespace eckit {
class PathName;
}

// ----------------------------------------------------------------------------

namespace atlas {
namespace option {

// ----------------------------------------------------------------------------

enum class FFT
{
    OFF    = 0,
    FFT992 = 1,
    FFTW   = 2,
};

// ----------------------------------------------------------------------------

class scalar_derivatives : public util::Config {
public:
    scalar_derivatives(bool);
};

// ----------------------------------------------------------------------------

class wind_EW_derivatives : public util::Config {
public:
    wind_EW_derivatives(bool);
};

// ----------------------------------------------------------------------------

class vorticity_divergence_fields : public util::Config {
public:
    vorticity_divergence_fields(bool);
};

// ----------------------------------------------------------------------------

class flt : public util::Config {
public:
    flt(bool);
};

// ----------------------------------------------------------------------------

class fft : public util::Config {
public:
    fft(FFT);
    fft(const std::string&);
};

class no_fft : public fft {
public:
    no_fft(): fft(FFT::OFF) {}
};

// ----------------------------------------------------------------------------

class split_latitudes : public util::Config {
public:
    split_latitudes(bool);
};

// ----------------------------------------------------------------------------

class write_legendre : public util::Config {
public:
    write_legendre(const eckit::PathName&);
};

// ----------------------------------------------------------------------------

class read_legendre : public util::Config {
public:
    read_legendre(const eckit::PathName&);
};

// ----------------------------------------------------------------------------

class write_fft : public util::Config {
public:
    write_fft(const eckit::PathName&);
};

// ----------------------------------------------------------------------------

class read_fft : public util::Config {
public:
    read_fft(const eckit::PathName&);
};

// ----------------------------------------------------------------------------

class nproma : public util::Config {
public:
    nproma(int);
};

// ----------------------------------------------------------------------------

class warning : public util::Config {
public:
    warning(int = 1);
};

// ----------------------------------------------------------------------------

}  // namespace option
}  // namespace atlas
