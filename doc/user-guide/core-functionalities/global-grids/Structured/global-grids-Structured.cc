#include "atlas/atlas.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Log.h"

using namespace atlas;
using atlas::grid::StructuredGrid;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    StructuredGrid grid( "O32" );

    Log::info() << "nlat   = " << grid.ny()  << std::endl;
    Log::info() << "nlon   = " << grid.nx(0) << std::endl;
    Log::info() << "npts   = " << grid.npts()  << std::endl;
    Log::info() << "lat    = " << grid.y(0)   << std::endl;
    Log::info() << "lon    = " << grid.x(1,0) << std::endl;
    Log::info() << "lonlat = " << grid.lonlat(1,0) << std::endl;

    atlas_finalize();
    return 0;
}



