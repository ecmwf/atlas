include(FetchContent)
set( HIC_ECBUILD_VERSION 3.13.1 )
message(STATUS "Downloading ecbuild version ${HIC_ECBUILD_VERSION} to ${CMAKE_BINARY_DIR}/ecbuild")
FetchContent_Populate(
    ecbuild
    URL            https://github.com/ecmwf/ecbuild/archive/refs/tags/${HIC_ECBUILD_VERSION}.tar.gz
    SOURCE_DIR     ${CMAKE_BINARY_DIR}/ecbuild
    QUIET
  )
find_package( ecbuild ${HIC_ECBUILD_VERSION} REQUIRED HINTS ${CMAKE_BINARY_DIR} )
