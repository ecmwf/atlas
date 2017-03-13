#include "atlas/library/atlas.h"
#include "atlas/runtime/Log.h"

int main(int argc, char** argv)
{
    atlas::atlas_init(argc, argv);
    atlas::Log::info() << "Hello world!" << std::endl;
    atlas::atlas_finalize();

    return 0;
}
