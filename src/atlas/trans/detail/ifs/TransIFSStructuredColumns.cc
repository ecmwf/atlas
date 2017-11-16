/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/trans/detail/ifs/TransIFSStructuredColumns.h"

#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/functionspace/Spectral.h"

namespace atlas {
namespace trans {

TransIFSStructuredColumns::TransIFSStructuredColumns(
    const functionspace::StructuredColumns& gp,
    const functionspace::Spectral& sp,
    const eckit::Configuration& config ) :
  TransIFS( gp.grid(), sp.truncation(), config ) {
}

TransIFSStructuredColumns::TransIFSStructuredColumns(
    const TransCache& cache,
    const functionspace::StructuredColumns& gp,
    const functionspace::Spectral& sp,
    const eckit::Configuration& config ) :
  TransIFS( cache, gp.grid(), sp.truncation(), config ) {
}


TransIFSStructuredColumns::~TransIFSStructuredColumns()
{
}

namespace {
static TransBuilderFunctionSpace< TransIFSStructuredColumns > builder("ifs(StructuredColumns,Spectral)");
}

} // namespace trans
} // namespace atlas
