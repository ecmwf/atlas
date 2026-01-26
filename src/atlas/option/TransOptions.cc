/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <map>

#include "eckit/filesystem/PathName.h"

#include "atlas/option/TransOptions.h"

// ----------------------------------------------------------------------------

namespace atlas {
namespace option {

vorticity_divergence_fields::vorticity_divergence_fields(bool v) {
    set("vorticity_divergence_fields", v);
}

wind_EW_derivatives::wind_EW_derivatives(bool v) {
    set("wind_EW_derivatives", v);
}

scalar_derivatives::scalar_derivatives(bool v) {
    set("scalar_derivatives", v);
}

flt::flt(bool flt) {
    set("flt", flt);
}

fft::fft(FFT fft) {
    static const std::map<FFT, std::string> FFT_to_string = {
        {FFT::OFF, "OFF"}, {FFT::FFT992, "FFT992"}, {FFT::FFTW, "FFTW"}, {FFT::pocketfft, "pocketfft"}};
    set("fft", FFT_to_string.at(fft));
}

fft::fft(const std::string& fft) {
    set("fft", fft);
}

split_y::split_y(bool split_y) {
    set("split_y", split_y);
}

write_legendre::write_legendre(const eckit::PathName& filepath) {
    set("write_legendre", filepath);
}

read_legendre::read_legendre(const eckit::PathName& filepath) {
    set("read_legendre", filepath);
}

write_fft::write_fft(const eckit::PathName& filepath) {
    set("write_fft", filepath);
}

read_fft::read_fft(const eckit::PathName& filepath) {
    set("read_fft", filepath);
}

nproma::nproma(int nproma) {
    set("nproma", nproma);
}

warning::warning(int warning) {
    set("warning", warning);
}

// ----------------------------------------------------------------------------

}  // namespace option
}  // namespace atlas
