#ifndef eckit_GribWrite_h
#define eckit_GribWrite_h

#include "grib_api.h"

#include "atlas/Mesh.hpp"
#include "atlas/grid/Field.h"

//------------------------------------------------------------------------------------------------------

namespace eckit {

//------------------------------------------------------------------------------------------------------

class GribWrite {
public:

    static void write( atlas::grid::FieldH& field, grib_handle* input_h );

    static void clone( atlas::grid::FieldH& field, const std::string& source, const std::string& output );

    static grib_handle* clone( atlas::grid::FieldH& field, grib_handle* source );

};

//---------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif

