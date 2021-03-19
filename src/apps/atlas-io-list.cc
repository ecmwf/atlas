/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>
#include <string>


#include "eckit/filesystem/PathName.h"

#include "atlas/io/Exceptions.h"
#include "atlas/io/RecordPrinter.h"
#include "atlas/io/print/Bytes.h"
#include "atlas/runtime/AtlasTool.h"

namespace atlas {


//----------------------------------------------------------------------------------------------------------------------

struct AtlasIOList : public atlas::AtlasTool {
    bool serial() override { return true; }
    int execute( const Args& args ) override;
    std::string briefDescription() override { return "Inspection of atlas-io files"; }
    std::string usage() override { return name() + " <file> [OPTION]... [--help,-h]"; }
    std::string longDescription() override {
        return "Inspection of atlas-io files\n"
               "\n"
               "       <file>: path to atlas-io file";
    }

    AtlasIOList( int argc, char** argv ) : AtlasTool( argc, argv ) {
        add_option( new SimpleOption<std::string>( "format", "Output format" ) );
        add_option( new SimpleOption<bool>( "version", "Print version of records" ) );
        add_option( new SimpleOption<bool>( "details", "Print detailed information" ) );
    }
};

//------------------------------------------------------------------------------------------------------

int AtlasIOList::execute( const Args& args ) {
    auto return_code = success();

    using namespace atlas;

    // User sanity checks
    if ( args.count() == 0 ) {
        Log::error() << "No file specified." << std::endl;
        help( std::cout );
        return failed();
    }

    // Configuration
    util::Config config;
    config.set( "format", args.getString( "format", "table" ) );
    config.set( "details", args.getBool( "details", false ) );

    // Loop over files
    for ( size_t f = 0; f < args.count(); ++f ) {
        eckit::PathName file( args( f ) );
        if ( !file.exists() ) {
            Log::error() << "File does not exist: " << file << std::endl;
            return failed();
        }
        auto filesize = size_t( file.size() );

        io::Session session;

        std::uint64_t pos = 0;
        try {
            while ( pos < filesize ) {
                auto uri    = io::Record::URI{file, pos};
                auto record = io::RecordPrinter{uri, config};

                std::stringstream out;
                out << "\n# " << uri.path << " [" << uri.offset << "]    "
                    << "{ size: " << atlas::io::Bytes{record.size()}.str( 0 ) << ",    version: " << record.version()
                    << ",    created: " << record.time() << " }";
                out << '\n' << ( config.getString( "format" ) == "table" ? "" : "---" ) << '\n';
                out << record << std::endl;

                std::cout << out.str();

                pos += record.size();
            }
        }
        catch ( const io::Exception& e ) {
            Log::error() << "    ATLAS-IO-ERROR: " << e.what() << std::endl;
            return_code = failed();
        }
    }
    return return_code;
}

}  // namespace atlas

//------------------------------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    atlas::AtlasIOList tool( argc, argv );
    return tool.start();
}
