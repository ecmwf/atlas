#ifndef atlas_Point3_h
#define atlas_Point3_h

#include "atlas/Parameters.hpp"

#include "eckit/container/sptree/SPPoint.h"

//-----------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

void latlon_to_3d( const double lat, const double lon, double* x, const double r, const double h );

void latlon_to_3d( const double lat, const double lon, double* x );

//-----------------------------------------------------------------------------

class Point3 : public eckit::SPPoint<3> {
public:

    Point3(): SPPoint<3>() {}

    double  operator[] (const size_t& i) const { assert(i<3); return x_[i]; }
    double& operator[] (const size_t& i)       { assert(i<3); return x_[i]; }

    Point3( double* p ): SPPoint<3>()
    {
        using namespace atlas;
        x_[XX] = p[XX];
        x_[YY] = p[YY];
        x_[ZZ] = p[ZZ];
    }

    Point3( const double lat, const double lon ): SPPoint<3>()
    {
        atlas::latlon_to_3d( lat, lon, x_ );
    }

    Point3( const double x, const double y, const double z ): SPPoint<3>(x,y)
    {
        x_[atlas::ZZ] = z;
    }

    template < typename T >
    void assign( const T& p )
    {
        using namespace atlas;
        x_[XX] = p[XX];
        x_[YY] = p[YY];
        x_[ZZ] = p[ZZ];
    }

    friend std::ostream& operator<<(std::ostream& s,const Point3& p);

};

//------------------------------------------------------------------------------------------------------

bool points_equal( const Point3& a, const Point3& b );

//---------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif

