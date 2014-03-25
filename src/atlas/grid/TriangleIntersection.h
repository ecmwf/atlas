#ifndef eckit_TriangleIntersection_h
#define eckit_TriangleIntersection_h

#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_DONT_ALIGN
#define EIGEN_DONT_VECTORIZE

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include "atlas/MeshGen.hpp"

#include "FloatCompare.h"

//-----------------------------------------------------------------------------

namespace eckit {

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

/// http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-9-ray-triangle-intersection/m-ller-trumbore-algorithm/

bool triag_intersection( const Triag& tg, const Ray& r, Isect& isect, const double slack = 2*std::numeric_limits<double>::epsilon() )
{
    Eigen::Vector3d edge1 = tg.v1 - tg.v0;
    Eigen::Vector3d edge2 = tg.v2 - tg.v0;
    Eigen::Vector3d pvec  = r.dir.cross(edge2);
    const double det = edge1.dot(pvec);
    if( std::abs(det) < slack ) return false;
    const double invDet = 1. / det;
    Eigen::Vector3d tvec = r.orig - tg.v0;
    isect.u = tvec.dot(pvec) * invDet;
    if( isect.u + slack < 0 || isect.u - slack > 1) return false;
    Eigen::Vector3d qvec = tvec.cross(edge1);
    isect.v = r.dir.dot(qvec) * invDet;
    if (isect.v + slack < 0 || isect.u + isect.v - slack > 1) return false;
    isect.t = edge2.dot(qvec) * invDet;
    return true;
}

//-----------------------------------------------------------------------------

} // namespace eckit

#endif

