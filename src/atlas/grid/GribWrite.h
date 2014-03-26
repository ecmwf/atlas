#ifndef eckit_GribWrite_h
#define eckit_GribWrite_h

#include "grib_api.h"

#include "atlas/Mesh.hpp"

//------------------------------------------------------------------------------------------------------

namespace eckit {

//------------------------------------------------------------------------------------------------------

class GribWrite {
public:

    static void write( atlas::Mesh& mesh, grib_handle* input_h );

    static void clone( atlas::Mesh& mesh, const std::string& source, const std::string& output );

    static grib_handle* clone( atlas::Mesh& mesh, grib_handle* source );

};

//---------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif

