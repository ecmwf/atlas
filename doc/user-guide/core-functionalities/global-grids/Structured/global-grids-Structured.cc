#include "atlas/atlas.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid/grids.h"
#include "eckit/config/Resource.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::grid::global;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    string gridID = Resource<string>( "--grid", string("N32") );

    SharedPtr<Structured> grid(Structured::create(gridID));

    Log::info() << "nlat   = " << grid->nlat()  << endl;
    Log::info() << "nlon   = " << grid->nlon(0) << endl;
    Log::info() << "npts   = " << grid->npts()  << endl;

    double lonlat[2];
    grid->lonlat(0, 1, lonlat);

    Log::info() << "lat    = " << grid->lat(0)   << endl;
    Log::info() << "lon    = " << grid->lon(0,1) << endl;
    Log::info() << "lonlat = " << lonlat[0] << "  "
                               << lonlat[1] << endl;

    atlas_finalize();
    return 0;
}



