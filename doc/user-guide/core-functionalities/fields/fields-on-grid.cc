#include "atlas/atlas.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid/grids.h"
#include "atlas/field/Field.h"
#include "eckit/config/Resource.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::array;
using namespace atlas::field;
using namespace atlas::grid::global;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    int jnode = 0;
    const double rpi = 2.0 * asin(1.0);
    const double deg2rad = rpi / 180.;
    const double zlatc = 0.0 * rpi;
    const double zlonc = 1.0 * rpi;
    const double zrad  = 2.0 * rpi / 9.0;
    double  zdist, zlon, zlat;

    string gridID;
    gridID = Resource<string>("--grid", string("N32"));
    SharedPtr<Structured> grid(Structured::create(gridID));
    int const nb_nodes = grid->npts();

    SharedPtr<Field> field_pressure(
      Field::create<double>("pressure", make_shape(nb_nodes)));

    ArrayView <double,1> pressure(*field_pressure);
    for (int jlat =0; jlat < grid->nlat(); ++jlat)
    {
        zlat = grid->lat(jlat);
        zlat = zlat * deg2rad;
        for (int jlon =0; jlon < grid->nlon(jlat); ++jlon)
        {
            zlon = grid->lon(jlat, jlon);
            zlon = zlon * deg2rad;
            zdist = 2.0 * sqrt((cos(zlat) * sin((zlon-zlonc)/2)) *
                              (cos(zlat) * sin((zlon-zlonc)/2)) +
                              sin((zlat-zlatc)/2) * sin((zlat-zlatc)/2));

            pressure(jnode) = 0.0;
            if (zdist < zrad)
            {
                pressure(jnode) = 0.5 * (1. + cos(rpi*zdist/zrad));
            }
            jnode = jnode+1;
        }
    }

    Log::info() << "==========================================" << endl;
    Log::info() << "memory field_pressure = "
                << field_pressure->bytes() * 1.e-9 << " GB"      << endl;
    Log::info() << "==========================================" << endl;

    atlas_finalize();

    return 0;
}
