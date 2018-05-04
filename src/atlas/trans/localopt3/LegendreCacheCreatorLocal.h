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

#include "atlas/trans/LegendreCacheCreator.h"
#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class LegendreCacheCreatorLocal : public trans::LegendreCacheCreatorImpl {
public:
    LegendreCacheCreatorLocal( const Grid&, int truncation, const eckit::Configuration& = util::NoConfig() );

    virtual ~LegendreCacheCreatorLocal();

    virtual std::string uid() const override;

    virtual void create(const std::string &path) const override;

private:
    Grid grid_;
    int truncation_;
    mutable std::string unique_identifier_;
};

// ------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
