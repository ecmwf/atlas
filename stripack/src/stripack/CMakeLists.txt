ecbuild_add_library(
  TARGET stripack
  SOURCES
   stripack.h
   stripack.cc
   stripack.f90
   stripack_interface.f90
  PUBLIC_INCLUDES
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
   $<INSTALL_INTERFACE:include>
   INSTALL_HEADERS ALL
   HEADER_DESTINATION  include/stripack 
)
