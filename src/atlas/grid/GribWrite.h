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

};

//---------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif

