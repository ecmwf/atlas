function( add_example target )
  set( options           "")
  set( single_value_args "")
  set( multi_value_args  "")
  cmake_parse_arguments(PARSE_ARGV 1 _PAR "${options}" "${single_value_args}" "${multi_value_args}")
  set( sources ${_PAR_UNPARSED_ARGUMENTS})
  ecbuild_add_executable( TARGET ${target} SOURCES ${sources} LIBS pluto NOINSTALL)
  pluto_declare_hic_source_files( ${sources} )
  pluto_target_linker_language_hic( ${target} )
endfunction()

add_example(pluto_example_is_on_device is_on_device.cc)
add_example(pluto_example_use_streams use_streams.cc)
add_example(pluto_example_use_unified_memory use_unified_memory.cc)
add_example(pluto_example_use_memory_pool use_memory_pool.cc)
add_example(pluto_example_use_scope use_scope.cc)

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module )

if (CMAKE_Fortran_COMPILER_ID MATCHES NVHPC AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 24.9 AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 25.5)
    # nvhpc 24.9 upto 25.5 have a compiler bug which fails to compile
    # pluto_example_mpl_allocator executable with errors like:
    #    NVFORTRAN-S-0155-Ambiguous interfaces for generic procedure allocate (pluto/examples/mpl_pinned_allocator.F90: 107)
    #    NVFORTRAN-S-0155-Ambiguous interfaces for generic procedure deallocate (pluto/examples/mpl_pinned_allocator.F90: 108)
    # Tested versions:
    #    working            : 22.11 24.5 24.7 25.5 25.7 25.9 25.11 26.1
    #    broken compilation : 24.9 24.11 25.1 25.3
    set(skip_pluto_example_mpl_allocator TRUE)
endif()
ecbuild_add_executable( TARGET pluto_example_mpl_allocator
                        SOURCES mpl_allocator.F90
                        LIBS pluto_f
                        NOINSTALL
                        CONDITION NOT skip_pluto_example_mpl_allocator )

