/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/option/Options.h"

// ----------------------------------------------------------------------------

namespace atlas {
namespace option {

vorticity_divergence_fields::vorticity_divergence_fields(bool v) { set("vorticity_divergence_fields",v); }

wind_EW_derivatives::wind_EW_derivatives(bool v) { set("wind_EW_derivatives",v); }

scalar_derivatives::scalar_derivatives(bool v) { set("scalar_derivatives",v); }

radius::radius(double _radius)
{
  set("radius",_radius);
}

radius::radius(const std::string &key)
{
  if( key == "Earth" ) {
    set("radius",util::Earth::radiusInMeters());
  } else {
    NOTIMP;
  }
}

type::type(const std::string &_type)
{
  set("type",_type);
}

halo::halo(size_t size)
{
  set("halo",size);
}

datatype::datatype(array::DataType::kind_t kind)
{
  set("datatype",kind);
}

datatype::datatype(const std::string &str )
{
  set("datatype",array::DataType::str_to_kind(str));
}

datatype::datatype(array::DataType dtype)
{
  set("datatype",dtype.kind());
}

name::name(const std::string &_name)
{
  set("name",_name);
}

global::global(size_t _owner)
{
  set("global",true);
  set("owner",_owner);
}

levels::levels(size_t _levels)
{
  set("levels",_levels);
}

variables::variables(size_t _variables)
{
  set("variables",_variables);
}

flt::flt(bool flt) { set("flt",flt); }

static const std::map<FFT,std::string> FFT_to_string = { {FFT992,"FFT992"},{FFTW,"FFTW"} };
//static const std::map<std::string,FFT> string_to_FFT = { {"FFT992",FFT992},{"FFTW",FFTW} };

fft::fft( FFT fft ) { set("fft",FFT_to_string.at(fft)); }
fft::fft( const std::string& fft ) { set("fft",fft); }
split_latitudes::split_latitudes(bool split_latitudes) { set("split_latitudes",split_latitudes); }
write_legendre::write_legendre( const eckit::PathName& filepath ) { set("write_legendre",filepath); }
read_legendre::read_legendre( const eckit::PathName& filepath ) { set("read_legendre",filepath); }


// ----------------------------------------------------------------------------

} // namespace option
} // namespace atlas
