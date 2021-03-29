/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "RecordPrinter.h"

#include <sstream>

#include "atlas/io/FileStream.h"
#include "atlas/io/print/JSONFormat.h"
#include "atlas/io/print/TableFormat.h"
#include "atlas/runtime/Exception.h"


namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

RecordPrinter::RecordPrinter( const eckit::PathName& path, const util::Config& config ) :
    RecordPrinter( path, 0, config ) {}

//---------------------------------------------------------------------------------------------------------------------

RecordPrinter::RecordPrinter( const eckit::PathName& path, const std::uint64_t offset, const util::Config& config ) :
    RecordPrinter( Record::URI{path, offset}, config ) {}

//---------------------------------------------------------------------------------------------------------------------

RecordPrinter::RecordPrinter( const Record::URI& ref, const util::Config& config ) :
    uri_( ref ), record_( Session::record( ref.path, ref.offset ) ) {
    if ( record_.empty() ) {
        auto in = InputFileStream( uri_.path );
        in.seek( uri_.offset );
        record_.read( in, true );
        ATLAS_ASSERT( not record_.empty() );
    }

    config.get( "format", options_.format );
    config.get( "details", options_.details );

    // Check if format is supported
    {
        std::vector<std::string> supported_formats{"json", "yaml", "table"};
        bool format_supported{false};
        for ( auto& supported_format : supported_formats ) {
            if ( options_.format == supported_format ) {
                format_supported = true;
                break;
            }
        }
        if ( not format_supported ) {
            std::stringstream s;
            s << "Format '" + options_.format + "' not supported. Supported formats:";
            for ( auto& supported_format : supported_formats ) {
                s << "\n  - " << supported_format;
            }
            throw_Exception( s.str(), Here() );
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

void RecordPrinter::print( std::ostream& out ) const {
    if ( options_.format == "json" ) {
        JSONFormat{uri_, util::Config( "details", options_.details )}.print( out );
    }
    else if ( options_.format == "yaml" ) {
        JSONFormat{uri_, util::Config( "details", options_.details )}.print( out );
    }
    else if ( options_.format == "table" ) {
        TableFormat{uri_, util::Config( "details", options_.details )}.print( out );
    }
    else {
        ATLAS_THROW_EXCEPTION( "Cannot print record: Unrecognized format " << options_.format << "." );
    }
}

//---------------------------------------------------------------------------------------------------------------------

std::ostream& operator<<( std::ostream& out, const RecordPrinter& info ) {
    info.print( out );
    return out;
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
