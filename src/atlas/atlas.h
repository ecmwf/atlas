#ifndef atlas_atlas_h
#define atlas_atlas_h

//#include "atlas_config.h"
//#include "atlas_defines.h"
//#include "atlas_version.h"

namespace atlas {

void atlas_init(int argc, char** argv);

extern "C" {
void atlas_initf(int argc, char** argv);
}

} // namespace atlas

#endif
