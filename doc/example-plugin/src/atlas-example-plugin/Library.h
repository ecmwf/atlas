/*
 * (C) Copyright 2023- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/Plugin.h"

namespace atlas {
namespace example_plugin {


class Library : public Plugin {
public:
    Library();
    static const Library& instance();

    std::string version() const override;
    std::string gitsha1( unsigned int count ) const override;
    void init() override;

};


}  // namespace example_plugin 
}  // namespace atlas
