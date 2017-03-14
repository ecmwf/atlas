#include "atlas/library/atlas.h"
#include "atlas/runtime/Log.h"

int main(int argc, char** argv)
{
    atlas::init(argc, argv);
    atlas::Log::info() << "Hello world!" << std::endl;
    atlas::finalise();

    return 0;
}
