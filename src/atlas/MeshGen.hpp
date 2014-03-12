// (C) Copyright 1996-2014 ECMWF.

#ifndef MeshGen_hpp
#define MeshGen_hpp

#include <cassert>
#include <string>

#include "atlas/Parameters.hpp"

namespace atlas {

class Mesh;

//------------------------------------------------------------------------------------------------------

/// Represents a 3D point
/// @note Do not assume memory to be initialised

struct Point3
{
    double x[3];

    double* data() { return x; }

    double operator()( const size_t& i ) const { assert( i < 3 ); return x[i]; }

    friend std::ostream& operator<<(std::ostream& s, const Point3& p )
    {
        s << '(' << p.x[XX] << "," << p.x[YY] << ","  << p.x[ZZ] << ')';
        return s;
    }
};

//------------------------------------------------------------------------------------------------------

void latlon_to_3d( const double lat, const double lon, double* x, const double r, const double h );
void latlon_to_3d( const double lat, const double lon, double* x );

//------------------------------------------------------------------------------------------------------

class MeshGen {

public:

    /// generate a mesh by triangulating the convex hull of the 3D points
    static Mesh* generate_from_points( const std::vector<Point3>& pts );

    /// generate regular spaced lat-long points
    static std::vector< Point3 >* generate_latlon_points( size_t nlats, size_t nlong );

    /// generates the cell centres en each cell
    /// @warning only for triangles ATM
    /// @warning this function must be checked for the correct INDEX translation to Fortran
    static void create_cell_centres( atlas::Mesh& mesh );
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // MeshGen_hpp
