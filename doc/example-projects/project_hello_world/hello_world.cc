#include "atlas/library.h"
#include "atlas/runtime/Log.h"

using atlas::Log;

int main( int argc, char* argv[] ) {
    atlas::initialize( argc, argv );
    Log::info() << "Hello world" << std::endl;
    Log::info() << "atlas version: " << atlas::Library::instance().version() << std::endl;
    atlas::finalize();
    return 0;
}
