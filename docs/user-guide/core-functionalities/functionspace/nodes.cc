#include "eckit/types/Types.h"
#include "eckit/memory/ScopedPtr.h"
#include "atlas/atlas.h"
#include "atlas/util/Debug.h"
#include "atlas/util/ArrayView.h"
#include "atlas/functionspace/Nodes.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/Mesh.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/Grid.h"
#include "atlas/Field.h"
#include "atlas/grids/ReducedGaussianGrid.h"


using namespace std;
using namespace atlas;
using namespace eckit;
using namespace atlas::grids;
using namespace atlas::meshgen;
using namespace atlas::functionspace;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);
    


    atlas_finalize();

    return 0;
}
