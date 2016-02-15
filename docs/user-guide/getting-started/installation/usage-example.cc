#include "atlas/atlas.h"

int main(int argc, char** argv)
{
    atlas::atlas_init(argc, argv);
    atlas::atlas_finalize();

    return 0;
}
