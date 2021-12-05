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

#include "atlas/trans/ifs/TransIFS.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
namespace functionspace {
class StructuredColumns;
class Spectral;
}  // namespace functionspace
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class TransIFSStructuredColumns : public trans::TransIFS {
public:
    TransIFSStructuredColumns(const functionspace::StructuredColumns&, const functionspace::Spectral&,
                              const eckit::Configuration& = util::Config());

    TransIFSStructuredColumns(const Cache&, const functionspace::StructuredColumns&, const functionspace::Spectral&,
                              const eckit::Configuration& = util::Config());

    virtual ~TransIFSStructuredColumns();
};

//-----------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
