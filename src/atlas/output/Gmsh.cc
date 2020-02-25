/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>

#include "atlas/output/Gmsh.h"
#include "atlas/output/detail/GmshImpl.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace output {

// -----------------------------------------------------------------------------

std::string GmshFileStream::parallelPathName( const eckit::PathName& path, int part ) {
    std::stringstream s;
    // s << path.dirName() << "/" << path.baseName(false) << "_p" << part <<
    // ".msh";
    s << path.asString() << ".p" << part;
    return s.str();
}

// -----------------------------------------------------------------------------

GmshFileStream::GmshFileStream( const eckit::PathName& file_path, const char* mode, int part ) : std::ofstream() {
    eckit::PathName par_path( file_path );
    std::ios_base::openmode omode = std::ios_base::out;
    if ( std::string( mode ) == "w" ) {
        omode = std::ios_base::out;
    }
    else if ( std::string( mode ) == "a" ) {
        omode = std::ios_base::app;
    }

    if ( part < 0 || mpi::size() == 1 ) {
        std::ofstream::open( file_path.localPath(), omode );
    }
    else {
        if ( mpi::rank() == 0 ) {
            eckit::PathName par_path( file_path );
            std::ofstream par_file( par_path.localPath(), std::ios_base::out );
            for ( idx_t p = 0; p < mpi::size(); ++p ) {
                par_file << "Merge \"" << parallelPathName( file_path, p ) << "\";" << std::endl;
            }
            par_file.close();
        }
        eckit::PathName path( parallelPathName( file_path, part ) );
        std::ofstream::open( path.localPath(), omode );
    }
}

//----------------------------------------------------------------------------------------------------------------------

Gmsh::Gmsh( const Output& output ) : Output( output ) {}

Gmsh::Gmsh( std::ostream& s ) : Output( new detail::GmshImpl( s ) ) {}

Gmsh::Gmsh( std::ostream& s, const eckit::Parametrisation& c ) : Output( new detail::GmshImpl( s, c ) ) {}

Gmsh::Gmsh( const eckit::PathName& p, const std::string& mode ) : Output( new detail::GmshImpl( p, mode ) ) {}

Gmsh::Gmsh( const eckit::PathName& p, const std::string& mode, const eckit::Parametrisation& c ) :
    Output( new detail::GmshImpl( p, mode, c ) ) {}

Gmsh::Gmsh( const eckit::PathName& p ) : Output( new detail::GmshImpl( p ) ) {}

Gmsh::Gmsh( const eckit::PathName& p, const eckit::Parametrisation& c ) : Output( new detail::GmshImpl( p, c ) ) {}

}  // namespace output
}  // namespace atlas
