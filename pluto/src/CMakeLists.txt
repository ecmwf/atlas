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
install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/pluto/pluto_config.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/pluto
)

ecbuild_add_library( TARGET pluto
  SOURCES
  #!${pluto_host_device_src}
    pluto/device/MemoryResource.cc
    pluto/device/allocator.h
    pluto/device/make_copy.h
    pluto/device/MemoryResource.h
    pluto/device/unique_ptr.h
    pluto/device/vector.h

    pluto/host/allocator.h
    pluto/host/make_copy.h
    pluto/host/MemoryResource.h
    pluto/host/unique_ptr.h
    pluto/host/vector.h

    pluto/memory_resource/compat/memory_resource
    pluto/memory_resource/yakl/YAKL_Gator.h
    pluto/memory_resource/yakl/YAKL_LinearAllocator.h
    pluto/memory_resource/DeviceMemoryResource.cc
    pluto/memory_resource/DeviceMemoryResource.h
    pluto/memory_resource/GatorMemoryResource.cc
    pluto/memory_resource/GatorMemoryResource.h
    pluto/memory_resource/ManagedMemoryResource.cc
    pluto/memory_resource/ManagedMemoryResource.h
    pluto/memory_resource/memory_resource.cc
    pluto/memory_resource/memory_resource.h
    pluto/memory_resource/MemoryPoolResource.cc
    pluto/memory_resource/MemoryPoolResource.h
    pluto/memory_resource/MemoryResourceAdaptor.h
    pluto/memory_resource/PinnedMemoryResource.cc
    pluto/memory_resource/PinnedMemoryResource.h
    pluto/memory_resource/TraceMemoryResource.cc
    pluto/memory_resource/TraceMemoryResource.h

    pluto/offload/copy.h
    pluto/offload/Event.cc
    pluto/offload/Event.h
    pluto/offload/memcpy.cc
    pluto/offload/memcpy.h
    pluto/offload/prefetch.cc
    pluto/offload/prefetch.h
    pluto/offload/Stream.cc
    pluto/offload/Stream.h
    pluto/offload/wait.cc
    pluto/offload/wait.h
    
    pluto/util/PointerInfo.cc
    pluto/util/PointerInfo.h
    
    pluto/util/Alignment.h
    pluto/util/PointerInfo.cc
    pluto/util/PointerInfo.h
    pluto/util/Registry.h
    pluto/util/Runtime.cc
    pluto/util/Runtime.h
    pluto/util/Trace.h

    pluto/pluto.h
    ${CMAKE_CURRENT_BINARY_DIR}/pluto/pluto_config.h

    INSTALL_HEADERS ALL
    HEADER_DESTINATION include
)

target_compile_features( pluto INTERFACE cxx_std_17 )
target_include_directories( pluto PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> )
target_link_libraries( pluto PUBLIC hic )
target_compile_options( pluto PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)


ecbuild_add_library( TARGET pluto_f 
		     SOURCES  pluto_f/pluto_module.F90
		              pluto_f/pluto_module.cc
		     PRIVATE_LIBS pluto
	           ) 
