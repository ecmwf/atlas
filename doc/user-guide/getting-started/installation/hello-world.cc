#include "atlas/library/Library.h"
#include "atlas/runtime/Log.h"

int main( int argc, char** argv ) {
    atlas::Library::instance().initialise( argc, argv );
    atlas::Log::info() << "Hello world!" << std::endl;
    atlas::Library::instance().finalise();

    return 0;
}
