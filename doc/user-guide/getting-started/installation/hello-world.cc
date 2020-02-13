#include "atlas/library.h"
#include "atlas/runtime/Log.h"

int main( int argc, char** argv ) {
    atlas::library::initialise( argc, argv );
    atlas::Log::info() << "Hello world!" << std::endl;
    atlas::library::finalise();

    return 0;
}
