/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/trans/detail/ifs/TransIFSNodeColumns.h"

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"

namespace atlas {
namespace trans {

TransIFSNodeColumns::TransIFSNodeColumns(
    const functionspace::NodeColumns& gp,
    const functionspace::Spectral& sp,
    const eckit::Configuration& ) :
  TransIFS( gp.mesh().grid(), sp.truncation() ) {
}


TransIFSNodeColumns::~TransIFSNodeColumns()
{
}


namespace {
static TransBuilderFunctionSpace< TransIFSNodeColumns > builder("ifs(NodeColumns,Spectral)");
}


} // namespace trans
} // namespace atlas
