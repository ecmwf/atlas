#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

using atlas::Log;
using atlas::StructuredGrid;

int main( int argc, char* argv[] ) {
    atlas::initialize( argc, argv );

    StructuredGrid grid( "O32" );

    Log::info() << "nx first  = " << grid.nx().front() << std::endl;
    Log::info() << "ny        = " << grid.ny() << std::endl;
    Log::info() << "npts      = " << grid.size() << std::endl;
    Log::info() << "xy        = " << grid.xy( 0, 0 ) << std::endl;
    Log::info() << "lonlat    = " << grid.lonlat( 0, 0 ) << std::endl;

    atlas::finalize();
    return 0;
}
