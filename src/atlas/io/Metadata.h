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

#include "atlas/io/Stream.h"
#include "atlas/io/detail/Checksum.h"
#include "atlas/io/detail/DataInfo.h"
#include "atlas/io/detail/Endian.h"
#include "atlas/io/detail/Link.h"
#include "atlas/io/detail/RecordInfo.h"
#include "atlas/io/detail/Type.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace io {

class Metadata;
class Stream;

//---------------------------------------------------------------------------------------------------------------------

size_t uncompressed_size( const atlas::io::Metadata& m );

//---------------------------------------------------------------------------------------------------------------------

class Metadata : public util::Config {
public:
    using util::Config::Config;

    Link link() const { return Link{getString( "link", "" )}; }

    Type type() const { return Type{getString( "type", "" )}; }

    void link( atlas::io::Metadata&& );

    std::string json() const;

    DataInfo data;
    RecordInfo record;
};

//---------------------------------------------------------------------------------------------------------------------

void write( const atlas::io::Metadata&, std::ostream& out );

void write( const atlas::io::Metadata&, atlas::io::Stream& out );

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
