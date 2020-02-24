#include "atlas/library.h"
#include "atlas/runtime/Log.h"

using atlas::Log;

int main( int argc, char* argv[] ) {
    atlas::library::initialise( argc, argv );
    Log::info() << "Hello world" << std::endl;
    Log::info() << "atlas version: " << atlas::Library::instance().version() << std::endl;
    atlas::library::finalise();
    return 0;
}
