#include "filesystem.h"

// POSIX dependencies
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "atlas/runtime/Exception.h"

namespace atlas::filesystem {

static void create_directory(const char path[]) {
    if (::mkdir(path, 0777) != 0) {
        ATLAS_THROW_EXCEPTION("Could not create directory: " << path);
    }
}

void create_directory(const std::string& path) {
    return create_directory(path.c_str());
}

static bool exists(const char path[]) {
    struct stat info;
    return stat( path, &info ) == 0;
}

bool exists(const std::string& path) {
    return exists(path.c_str());
}

static void remove_all(const char path[]) {
    size_t path_len;
    char *full_path;
    DIR *dir;
    struct stat stat_path, stat_entry;
    struct dirent *entry;

    // stat for the path
    ::stat(path, &stat_path);

    // if path does not exists or is not dir, throw
    if (S_ISDIR(stat_path.st_mode) == 0) {
        ATLAS_THROW_EXCEPTION("Is not a directory: " << path);
    }

    // if not possible to read the directory for this user
    if ((dir = opendir(path)) == nullptr) {
        ATLAS_THROW_EXCEPTION("Cannot open directory: " << path);
    }

    // the length of the path
    path_len = ::strlen(path);

    // iteration through entries in the directory
    while ((entry = readdir(dir)) != nullptr) {

        // skip entries "." and ".."
        if (!::strcmp(entry->d_name, ".") || !::strcmp(entry->d_name, "..")) {
            continue;
        }

        // determinate a full path of an entry
        full_path = (char*) ::calloc(path_len + 1 + ::strlen(entry->d_name) + 1, sizeof(char));
        ::strcpy(full_path, path);
        ::strcat(full_path, "/");
        ::strcat(full_path, entry->d_name);

        // stat for the entry
        ::stat(full_path, &stat_entry);

        // recursively remove a nested directory
        if (S_ISDIR(stat_entry.st_mode) != 0) {
            remove_all(full_path);
            ::free(full_path);
            continue;
        }

        // remove a file object
        if (::unlink(full_path) != 0) {
            ATLAS_THROW_EXCEPTION("Can't remove a file: " << full_path);
        }
        ::free(full_path);
    }

    // remove the devastated directory and close the object of it
    if (::rmdir(path) != 0) {
        ATLAS_THROW_EXCEPTION("Can't remove a directory: " << path);
    }
    ::closedir(dir);
}

void remove_all(const std::string& path) {
    remove_all(path.c_str());
}

}
