#include "atlas/grid/PointIndex3.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

PointIndex3* create_cell_centre_index( atlas::Mesh& mesh )
{
    atlas::FunctionSpace& triags = mesh.function_space( "triags" );
    atlas::FieldT<double>& triags_centres = triags.field<double>( "centre" );

    const size_t npts = triags.bounds()[1];

    std::vector<typename PointIndex3::Value> p;
    p.reserve(npts);

    for( size_t ip = 0; ip < npts; ++ip )
    {
        p.push_back( typename PointIndex3::Value(
                         typename PointIndex3::Point( triags_centres(atlas::XX,ip),
                                                      triags_centres(atlas::YY,ip),
                                                      triags_centres(atlas::ZZ,ip) ), ip ) );
    }

    PointIndex3* tree = new PointIndex3();

    tree->build(p.begin(), p.end());

    return tree;
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

