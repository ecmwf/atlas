/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/trans/ifs/TransIFSStructuredColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/trans/detail/TransFactory.h"

namespace atlas {
namespace trans {

TransIFSStructuredColumns::TransIFSStructuredColumns( const functionspace::StructuredColumns& gp,
                                                      const functionspace::Spectral& sp,
                                                      const eckit::Configuration& config ) :
    TransIFSStructuredColumns( Cache(), gp, sp, config ) {}

TransIFSStructuredColumns::TransIFSStructuredColumns( const Cache& cache, const functionspace::StructuredColumns& gp,
                                                      const functionspace::Spectral& sp,
                                                      const eckit::Configuration& config ) :
    TransIFS( cache, gp.grid(), sp.truncation(), config ) {
    assertCompatibleDistributions( gp, sp );
    spectral_ = sp;
}

TransIFSStructuredColumns::~TransIFSStructuredColumns() {}

namespace {
static TransBuilderFunctionSpace<TransIFSStructuredColumns> builder( "ifs(StructuredColumns,Spectral)", "ifs" );
}

}  // namespace trans
}  // namespace atlas
