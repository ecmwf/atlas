list( APPEND atlas_testing_src
  UnitTest.hpp
  Macros.hpp
  dummy.cpp
)

# add_library(atlas_testing STATIC ${atlas_testing_src})
# set_target_properties(atlas_testing PROPERTIES LINKER_LANGUAGE CXX)

ecbuild_add_library(TARGET atlas_testing
          AUTO_VERSION
          TYPE STATIC
          LINKER_LANGUAGE CXX
          INSTALL_HEADERS ALL
          SOURCES
            ${atlas_testing_src}
          HEADER_DESTINATION
            ${INSTALL_INCLUDE_DIR}/atlas/testing
          )
