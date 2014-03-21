// (C) Copyright 1996-2014 ECMWF.

#ifndef MeshGen_hpp
#define MeshGen_hpp

#include <cassert>
#include <string>
#include <ostream>

#include "atlas/Parameters.hpp"

namespace atlas {

class Mesh;

//------------------------------------------------------------------------------------------------------

/// Represents a 3D point
/// @note Do not assume memory to be initialised
/// @todo remove this class, algorithm should work with existing coordinates on mesh

struct Point3
{
    Point3() {}

    Point3( double* p )
    {
        x_[XX] = p[XX];
        x_[YY] = p[YY];
        x_[ZZ] = p[ZZ];
    }

    double* data() { return x_; }

    double operator()( const size_t& i ) const { assert( i < 3 ); return x_[i]; }

    friend std::ostream& operator<<(std::ostream& s, const Point3& p )
    {
        s << '(' << p.x_[XX] << "," << p.x_[YY] << ","  << p.x_[ZZ] << ')';
        return s;
    }

private:
    double x_[3];
};

//------------------------------------------------------------------------------------------------------

void latlon_to_3d( const double lat, const double lon, double* x, const double r, const double h );
void latlon_to_3d( const double lat, const double lon, double* x );

//------------------------------------------------------------------------------------------------------

class MeshGen {

public:

    /// generate a mesh by triangulating the convex hull of the 3D points
    static Mesh* generate_from_points( const std::vector<Point3>& pts );

    /// generate regular spaced lat-long points (does not include extremes)
    static std::vector< Point3 >* generate_latlon_points( size_t nlats, size_t nlong );

    /// generate regular lat-long grid points (includes extremes -90,90 and 0,360)
    static std::vector< Point3 >* generate_latlon_grid( size_t nlats, size_t nlong );

    /// generates the cell centres en each cell
    /// @warning only for triangles ATM
    /// @warning this function must be checked for the correct INDEX translation to Fortran
    static void create_cell_centres( atlas::Mesh& mesh );
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // MeshGen_hpp
