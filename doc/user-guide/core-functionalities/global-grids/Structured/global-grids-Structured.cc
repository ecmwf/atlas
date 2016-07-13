#include "atlas/atlas.h"
#include "atlas/grid/Structured.h"
#include "atlas/runtime/Log.h"

using atlas::atlas_init;
using atlas::atlas_finalize;
using atlas::grid::Structured;
using atlas::Log;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    Structured::Ptr grid(Structured::create( "O32" ));

    Log::info() << "nlat   = " << grid->nlat()  << std::endl;
    Log::info() << "nlon   = " << grid->nlon(0) << std::endl;
    Log::info() << "npts   = " << grid->npts()  << std::endl;

    double lonlat[2];
    grid->lonlat(0, 1, lonlat);

    Log::info() << "lat    = " << grid->lat(0)   << std::endl;
    Log::info() << "lon    = " << grid->lon(0,1) << std::endl;
    Log::info() << "lonlat = " << lonlat[0] << "  "
                               << lonlat[1] << std::endl;

    atlas_finalize();
    return 0;
}



