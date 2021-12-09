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

#include "atlas/grid/Grid.h"
#include "atlas/trans/LegendreCacheCreator.h"
#include "atlas/util/Config.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class LegendreCacheCreatorLocal : public trans::LegendreCacheCreatorImpl {
public:
    LegendreCacheCreatorLocal(const Grid&, int truncation, const eckit::Configuration& = util::NoConfig());

    virtual ~LegendreCacheCreatorLocal() override;

    virtual bool supported() const override;

    virtual std::string uid() const override;

    virtual void create(const std::string& path) const override;

    virtual Cache create() const override;

    size_t estimate() const override;

private:
    const Grid grid_;
    const int truncation_;
    const util::Config config_;
    mutable std::string unique_identifier_;
};

// ------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
