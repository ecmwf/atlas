#include <cstring>
#include <fstream>
#include <sstream>

#include "eckit/filesystem/PathName.h"
#include "eckit/runtime/Main.h"

#include "atlas/runtime/Exception.h"
#include "atlas_f/internals/atlas_read_file.h"

namespace atlas {

void read_file( const eckit::PathName& p, std::ostream& out ) {
    if ( p.exists() ) {
        std::ifstream in;
        in.open( p.asString().c_str() );
        if ( !in ) { throw_CantOpenFile( p.asString(), Here() ); }
        else {
            out << in.rdbuf();
            in.close();
        }
    }
}

}  // namespace atlas

extern "C" {

//-----------------------------------------------------------------------------

int atlas__read_file( const char* path, char*& content, int& size ) {
    // eckit::FileReadPolicy p =
    // eckit::Main::instance().behavior().fileReadPolicy();

    // std::stringstream ss;

    // if( read( p, path, ss ) )
    // {
    //   std::string s = ss.str();
    //   size = s.size()+1;
    //   content = new char[size];
    //   strcpy(content,s.c_str());
    //   return true;
    // }

    std::stringstream ss;
    atlas::read_file( path, ss );
    std::string s = ss.str();
    size          = s.size() + 1;
    content       = new char[size];
    strcpy( content, s.c_str() );
    return true;
}

//-----------------------------------------------------------------------------
}
