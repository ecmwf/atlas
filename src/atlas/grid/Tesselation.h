// (C) Copyright 1996-2014 ECMWF.

#ifndef atlas_grid_Tesselation_hpp
#define atlas_grid_Tesselation_hpp

#include "atlas/Mesh.hpp"

namespace atlas {

//------------------------------------------------------------------------------------------------------

class Tesselation {

public:

    /// generate a mesh by triangulating the convex hull of the 3D points
    static void tesselate( atlas::Mesh& mesh );

    /// generate regular spaced lat-long points (does not include extremes)
    static void generate_latlon_points( atlas::Mesh& mesh, const size_t& nlats, const size_t& nlong );

    /// generate regular lat-long grid points (includes extremes -90,90 and 0,360)
    static void generate_latlon_grid( atlas::Mesh& mesh, const size_t& nlats, const size_t& nlong );

    /// generates the cell centres en each cell
    /// @warning only for triangles ATM
    /// @warning this function must be checked for the correct INDEX translation to Fortran
    static void create_cell_centres( atlas::Mesh& mesh );

    /// Ensures or creates the necessary nodal structure in the mesh object
    static void create_mesh_structure( atlas::Mesh& mesh, const size_t nb_nodes );

};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
