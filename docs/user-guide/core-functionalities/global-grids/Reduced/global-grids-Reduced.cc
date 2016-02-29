#include "atlas/atlas.h"
#include "atlas/grid/grids.h"
#include "eckit/config/Resource.h"

using namespace std;
using namespace atlas;
using namespace atlas::grids;

int main(int argc, char *argv[])
{    
    atlas_init(argc, argv);

    string gridID;
    gridID = eckit::Resource<string>("--grid", string("N32"));

    ReducedGrid::Ptr reducedGrid(ReducedGrid::create(gridID));

    cout << "nlat   = " << reducedGrid->nlat()  << endl;
    cout << "nlon   = " << reducedGrid->nlon(0) << endl;
    cout << "npts   = " << reducedGrid->npts()  << endl;

    double crdLat, crdLon, crdLonLat[2];
    crdLat = reducedGrid->lat(0);
    crdLon = reducedGrid->lon(0, 1);
    reducedGrid->lonlat(0, 1, crdLonLat);

    cout << "lat    = " << crdLat << endl;
    cout << "lon    = " << crdLon << endl;
    cout << "lonlat = " << crdLonLat[0] << "  "
                        << crdLonLat[1] << endl;

    atlas_finalize();
    return 0;
}



