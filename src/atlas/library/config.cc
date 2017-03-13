#include <algorithm>
#include <string>

#include "atlas/library/config.h"
#include "atlas/library/git_sha1.h"

namespace atlas {

const char* version() {
    return ATLAS_VERSION;
}

int version_int() {
    return 10000*ATLAS_MAJOR_VERSION + 100*ATLAS_MINOR_VERSION + 1*ATLAS_PATCH_VERSION;
}

const char* git_sha1(unsigned int chars) {
    static std::string sha1(ATLAS_GIT_SHA1);
    if(sha1.empty()) {
        return "not available";
    }
    sha1 = sha1.substr(0, std::min(chars,40u));
    return sha1.c_str();
}

}
