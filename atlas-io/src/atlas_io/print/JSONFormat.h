/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>
#include <map>
#include <string>

#include "eckit/config/Configuration.h"

#include "atlas_io/Metadata.h"
#include "atlas_io/Record.h"

namespace atlas {
namespace io {


class JSONFormat {
public:
    JSONFormat(const Record::URI& record, const eckit::Configuration&);

    void print(std::ostream&) const;

private:
    const Record record_;
    std::map<std::string, Metadata> items_;

    bool print_details_{false};
};


}  // namespace io
}  // namespace atlas
