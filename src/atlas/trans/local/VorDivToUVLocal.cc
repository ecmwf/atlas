/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/trans/local/VorDivToUVLocal.h"

#include "atlas/functionspace/Spectral.h"
#include "atlas/runtime/Log.h"

using atlas::functionspace::Spectral;
using atlas::FunctionSpace;

namespace atlas {
namespace trans {

namespace {
static VorDivToUVBuilder<VorDivToUVLocal> builder("local");
}

void VorDivToUVLocal::execute( const int nb_coeff, const int nb_fields,
                               const double vorticity[], const double divergence[],
                               double U[], double V[],
                               const eckit::Configuration& config ) const
{
  NOTIMP;
}


VorDivToUVLocal::VorDivToUVLocal( const int truncation, const eckit::Configuration& config ) :
  truncation_(truncation) {
}

VorDivToUVLocal::VorDivToUVLocal( const FunctionSpace& fs, const eckit::Configuration& config ) :
  truncation_( Spectral(fs).truncation() ) {
}


VorDivToUVLocal::~VorDivToUVLocal()
{
}

} // namespace trans
} // namespace atlas
