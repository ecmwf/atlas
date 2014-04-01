#ifndef atlas_grid_TriangleIntersection_h
#define atlas_grid_TriangleIntersection_h

#include "eckit/maths/Eigen.h"

#include "eckit/types/FloatCompare.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

/// Intersection data structure

struct Isect
{
    double u;
    double v;
    double t;

    double w() const { return 1.0 - u - v; }

    friend std::ostream& operator<<(std::ostream& s,const Isect& p)
    {
        s << '(' << p.u << "," << p.v << ","  << p.w() << ","  << p.t << ')';
        return s;
    }

};

/// Ray trace data structure

struct Ray
{
    Eigen::Vector3d orig;
    Eigen::Vector3d dir;

    /// initializes ray with origin in point and direction to (0,0,0)
    Ray( double* p )
    {
        orig = Eigen::Vector3d::Map(p);
        dir  = - orig;
    }

    Ray( double* o, double* d )
    {
        orig = Eigen::Vector3d::Map(o);
        dir  = Eigen::Vector3d::Map(d);
    }

    Eigen::Vector3d operator()( double t ) const { return orig + t*dir; }
};

/// triangle structure

struct Triag
{
    Triag( double* x0, double* x1, double* x2 )
    {
        v0 = Eigen::Vector3d::Map(x0);
        v1 = Eigen::Vector3d::Map(x1);
        v2 = Eigen::Vector3d::Map(x2);
    }

    Eigen::Vector3d v0;
    Eigen::Vector3d v1;
    Eigen::Vector3d v2;
};

//------------------------------------------------------------------------------------------------------

/// http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-9-ray-triangle-intersection/m-ller-trumbore-algorithm

bool triag_intersection( const Triag& tg, const Ray& r, Isect& isect, const double slack = 2*std::numeric_limits<double>::epsilon() );

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif

