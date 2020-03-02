#include "atlas/library.h"
#include "atlas/runtime/Log.h"

int main( int argc, char** argv ) {
    atlas::initialize( argc, argv );
    atlas::Log::info() << "Hello world!" << std::endl;
    atlas::finalize();

    return 0;
}
