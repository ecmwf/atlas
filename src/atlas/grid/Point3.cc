#include "atlas/grid/Point3.h"
#include "atlas/grid/FloatCompare.h"

using namespace eckit;

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& s,const Point3& p)
{
    s << '(' << p.x_[atlas::XX] << ","
             << p.x_[atlas::YY] << ","
             << p.x_[atlas::ZZ] << ')';
    return s;
}

//------------------------------------------------------------------------------------------------------

const double radius_earth = 6367.47; // from ECMWF model ...

void latlon_to_3d(const double lat, const double lon, double* x, const double r, const double h )
{
    // See http://en.wikipedia.org/wiki/Geodetic_system#From_geodetic_to_ECEF

    double& X = x[XX];
    double& Y = x[YY];
    double& Z = x[ZZ];

    const double a = r;   // 6378137.0 ;       // WGS84 semi-major axis
    const double e2 = 0;  // ignored -- 6.69437999014E-3; // WGS84 first numerical eccentricity squared

    const double phi = lat / 180.0 * M_PI;
    const double lambda = lon / 180.0 * M_PI;

    const double cos_phi = cos(phi);
    const double sin_phi = sin(phi);
    const double cos_lambda = cos(lambda);
    const double sin_lambda = sin(lambda);

    const double N_phi = a/sqrt(1-e2*sin_phi*sin_phi);

    X = (N_phi + h) * cos_phi * cos_lambda;
    Y = (N_phi + h) * cos_phi * sin_lambda;
    Z = (N_phi * (1-e2) + h) * sin_phi;
}

void latlon_to_3d( const double lat, const double lon, double* x )
{
    latlon_to_3d( lat, lon, x, radius_earth, 0. );
}

//------------------------------------------------------------------------------------------------------

bool points_equal(const Point3 &a, const Point3 &b)
{
    return FloatCompare::is_equal( Point3::distance2(a,b), 0.0 );
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

