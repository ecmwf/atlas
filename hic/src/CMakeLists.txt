# (C) Copyright 2024- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set(HIC_BACKEND_CUDA  0)
set(HIC_BACKEND_HIP   0)
set(HIC_BACKEND_DUMMY 0)
if( HAVE_CUDA )
  set(HIC_BACKEND_CUDA  1)
elseif( HAVE_HIP )
  set(HIC_BACKEND_HIP   1)
else()
  set(HIC_BACKEND_DUMMY 1)
endif()

configure_file(${PROJECT_SOURCE_DIR}/src/hic/hic_config.h.in ${PROJECT_BINARY_DIR}/src/hic/hic_config.h @ONLY )

ecbuild_add_library( TARGET hic
  TYPE INTERFACE
  PUBLIC_INCLUDES
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

install( FILES ${PROJECT_BINARY_DIR}/src/hic/hic_config.h DESTINATION include/hic )
install( FILES hic/hic.h                                  DESTINATION include/hic )
install( FILES hic/hic_runtime.h                          DESTINATION include/hic )
install( FILES hic/hic_dummy/hic_dummy_runtime.h          DESTINATION include/hic/hic_dummy )

if( HAVE_CUDA )
  target_link_libraries( hic INTERFACE CUDA::cudart )
elseif( HAVE_HIP )
  target_link_libraries( hic INTERFACE hip::host )
endif()
