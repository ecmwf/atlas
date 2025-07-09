# (C) Copyright 2013- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


#pluto_host_device( pluto_host_device_src SOURCE pluto/device/MemoryResource.cc)

set(PLUTO_HAVE_CUDA ${HAVE_CUDA} )
set(PLUTO_HAVE_HIP  ${HAVE_HIP} )
set(PLUTO_HAVE_HIC  0)

if( HAVE_CUDA OR HAVE_HIP )
  set(PLUTO_HAVE_HIC 1)
endif()

configure_file(pluto/pluto_config.h.in pluto/pluto_config.h @ONLY )
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/pluto/pluto_config.h
  DESTINATION
    ${INSTALL_INCLUDE_DIR}/pluto
)

if (NOT HAVE_PMR)
  list( APPEND pmr_compat_src
      pluto/memory_resource/compat/memory_resource
      pluto/memory_resource/compat/memory_resource.cc
  )
endif()

ecbuild_add_library( TARGET pluto
  SOURCES
    pluto/pluto.h
    ${CMAKE_CURRENT_BINARY_DIR}/pluto/pluto_config.h

    pluto/device/allocator.h
    pluto/device/MemoryResource.cc
    pluto/device/MemoryResource.h

    pluto/host/allocator.h
    pluto/host/MemoryResource.cc
    pluto/host/MemoryResource.h

    pluto/memory.cc
    pluto/memory.h
    pluto/memory_resource.h
    pluto/memory_resource/AsyncMemoryResourceAdaptor.cc
    pluto/memory_resource/AsyncMemoryResourceAdaptor.h
    pluto/memory_resource/DeviceMemoryResource.cc
    pluto/memory_resource/DeviceMemoryResource.h
    pluto/memory_resource/HostMemoryResource.cc
    pluto/memory_resource/HostMemoryResource.h
    pluto/memory_resource/ManagedMemoryResource.cc
    pluto/memory_resource/ManagedMemoryResource.h
    pluto/memory_resource/memory_resource.cc
    pluto/memory_resource/MemoryPoolResource.cc
    pluto/memory_resource/MemoryPoolResource.h
    pluto/memory_resource/MemoryResourceAdaptor.h
    pluto/memory_resource/PinnedMemoryResource.cc
    pluto/memory_resource/PinnedMemoryResource.h
    pluto/memory_resource/TraceMemoryResource.cc
    pluto/memory_resource/TraceMemoryResource.h
    pluto/memory_resource/detail/GatorMemoryResource.cc
    pluto/memory_resource/detail/GatorMemoryResource.h
    pluto/memory_resource/detail/yakl/YAKL_Gator.h
    pluto/memory_resource/detail/yakl/YAKL_LinearAllocator.h

    pluto/alignment.h
    pluto/copy.h
    pluto/event.cc
    pluto/event.h
    pluto/memcpy.cc
    pluto/memcpy.h
    pluto/pointer_info.cc
    pluto/pointer_info.h
    pluto/prefetch.cc
    pluto/prefetch.h
    pluto/runtime.cc
    pluto/runtime.h
    pluto/scope.cc
    pluto/scope.h
    pluto/stream.cc
    pluto/stream.h
    pluto/trace.cc
    pluto/trace.h
    pluto/wait.cc
    pluto/wait.h
    
    pluto/detail/Registry.h

    ${pmr_compat_src}

  INSTALL_HEADERS ALL
  HEADER_DESTINATION include
  PUBLIC_LIBS hic
  PUBLIC_INCLUDES
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

target_compile_features( pluto INTERFACE cxx_std_17 )
target_compile_options( pluto PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)

if (HAVE_FORTRAN)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/pluto)
  ecbuild_add_library( TARGET pluto_f
    SOURCES
      pluto_f/pluto_module.F90
      pluto_f/pluto_module.cc
    PRIVATE_LIBS pluto
    PUBLIC_INCLUDES $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
                    $<INSTALL_INTERFACE:module/pluto>
  )
  if( ECBUILD_INSTALL_FORTRAN_MODULES )
    install( DIRECTORY ${CMAKE_BINARY_DIR}/module/pluto/
             DESTINATION module/pluto
             COMPONENT modules )
  endif()
endif()

