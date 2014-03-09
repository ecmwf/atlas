// (C) Copyright 1996-2014 ECMWF.

#ifndef MeshGen_hpp
#define MeshGen_hpp

#include <string>

namespace atlas {

class Mesh;

//------------------------------------------------------------------------------------------------------

/// Represents a 3D point
/// @note Do not assume memory to be initialised

struct Point3
{
    double x[3];
};

//------------------------------------------------------------------------------------------------------

void latlon_to_3d( const double lat, const double lon, double* x, const double r, const double h );
void latlon_to_3d( const double lat, const double lon, double* x );

//------------------------------------------------------------------------------------------------------

class MeshGen {

public:

    /// generate a mesh by triangulating the convex hull of the 3D points
    static Mesh* generate_from_ll_points( const std::vector<Point3>& pts );

};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // MeshGen_hpp
