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

#include <cstdint>
#include <iosfwd>
#include <string>

#include "eckit/filesystem/PathName.h"

#include "atlas/io/Record.h"
#include "atlas/io/Session.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class RecordPrinter {
public:
    RecordPrinter( const Record::URI&, const util::Config& = util::NoConfig() );

    RecordPrinter( const eckit::PathName&, const util::Config& = util::NoConfig() );

    RecordPrinter( const eckit::PathName&, std::uint64_t offset, const util::Config& = util::NoConfig() );

    Record record() const { return record_; }

    size_t size() const { return record_.size(); }

    Version version() const { return record_.version(); }

    Time time() const { return record_.time(); }

    void print( std::ostream& out ) const;

    friend std::ostream& operator<<( std::ostream&, const RecordPrinter& );

private:
    Session session_;

    Record::URI uri_;

    struct {
        std::string format{"table"};
        bool details{false};
    } options_;

    Record record_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
