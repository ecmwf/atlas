
ecbuild_add_executable( TARGET pluto-benchmark 
    SOURCES run_benchmark.cc
    LIBS    pluto
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} )
target_host_device_sources( pluto-benchmark PRIVATE benchmark.cc)
