#pragma once

#include <string>

namespace atlas::filesystem {

// See c++ std::filesystem for API
void create_directory(const std::string& path);
bool exists(const std::string& path);
void remove_all(const std::string& path);

}
