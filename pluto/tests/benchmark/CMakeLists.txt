
ecbuild_add_executable( TARGET pluto-benchmark 
    SOURCES run_benchmark.cc benchmark.cc
    LIBS    pluto
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} )
pluto_declare_hic_source_files(benchmark.cc)
