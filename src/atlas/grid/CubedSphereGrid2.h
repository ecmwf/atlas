#pragma once

#include "atlas/grid/Grid.h"

namespace atlas {

class CubedSphereGrid2 : public atlas::Grid {
public:
    CubedSphereGrid2(idx_t resolution);
    void printCubedSphere2() const ; // Temporary
};

}  // namespace atlas
