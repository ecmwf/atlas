# (C) Copyright 2025 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


ecbuild_add_executable(
  TARGET  atlas-interpolations
  SOURCES atlas-interpolations.cc
          ScripIO.h
          ScripIO.cc
  LIBS    atlas
)

# Overcome Cray linking problem where the C++ library is not linked with the C library
list( APPEND NetCDF_CXX_EXTRA_LIBRARIES NetCDF::NetCDF_C )

ecbuild_add_option( FEATURE NETCDF
                    DESCRIPTION "Compile support the SCRIP format"
                    REQUIRED_PACKAGES "NetCDF COMPONENTS C CXX" )

target_compile_definitions( atlas-interpolations PRIVATE ATLAS_HAVE_NETCDF=${HAVE_NETCDF} )
if( HAVE_NETCDF )
  target_link_libraries( atlas-interpolations NetCDF::NetCDF_CXX )
endif()

function( atlas_test_app_atlas_interpolations case nprocs )

    string(REPLACE "/" ";" case ${case} )
    list( GET case 0 sgrid )
    list( GET case 1 tgrid )
    list( GET case 2 itype )
    unset( checksums )
    unset( case_tests )

    if( NOT DEFINED MPI_SLOTS )
      set( MPI_SLOTS 9999 )
    endif()

    foreach( n ${nprocs} )
      if ( (n EQUAL 0) OR eckit_HAVE_MPI )
        if( MPI_SLOTS GREATER_EQUAL ${n} )
            foreach( t 1 )
                set( case_name atlas_test_app_atlas_interpolations_${sgrid}_to_${tgrid}_with_${itype}_n${n}_t${t} )
                set( test_name ${case_name}_write )
                ecbuild_add_test( TARGET ${test_name}
                    COMMAND atlas-interpolations
                    ARGS    --display-name=${test_name}
                            --interpolate
                            --s.grid=${sgrid}
                            --t.grid=${tgrid}
                            --i.type=${itype}
                            --test-matrix
                            --output-matrix
                            --matrix.name=${case_name}
                            --matrix.format=eckit
                            --output-gmsh
                            --gmsh.coordinates=lonlat
                            --test-interpolator
                            --i.knn=4 # just in case itype==knn
                            --output-checksum
                            --checksum.file=${test_name}.checksum
                    MPI ${n}
                    OMP ${t}
                    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                )

                set( test_name ${case_name}_read )
                ecbuild_add_test( TARGET ${test_name}
                    COMMAND atlas-interpolations
                    ARGS    --display-name=${test_name}
                            --interpolate
                            --s.grid=${sgrid}
                            --t.grid=${tgrid}
                            --i.type=${itype}
                            --test-matrix
                            --read-matrix
                            --matrix.name=${case_name}
                            --matrix.format=eckit
                            --test-interpolator
                            --output-gmsh
                            --i.knn=4 # just in case itype==knn
                            --output-checksum
                            --checksum.file=${test_name}.checksum
                    MPI ${n}
                    OMP ${t}
                    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                )
                set_tests_properties( ${case_name}_read PROPERTIES DEPENDS ${case_name}_write )

                # list( APPEND checksums ${test_name}.checksum )
                list( APPEND case_tests ${test_name} )
            endforeach()
        endif()
      endif()
    endforeach()
    # This kind of checksum test can be added, inspired from in tests/acceptance_tests/atlas_atest_mgrids
    # add_test(
    #     NAME atlas_test_app_atlas_interpolations_${sgrid}_to_${tgrid}_with_${itype}_BitReproducible
    #     COMMAND ${CMAKE_COMMAND}
    #     "-DFILES=${checksums}" -P ${CMAKE_CURRENT_SOURCE_DIR}/compare.cmake
    # )
    # foreach( dep ${case_tests} )
    #     set_tests_properties( atlas_test_app_atlas_interpolations_${sgrid}_to_${tgrid}__with_${itype}_BitReproducible PROPERTIES DEPENDS ${dep} )
    # endforeach()
endfunction()

if (HAVE_TESTS)
  unset(cases)
  list( APPEND cases
    O16/O32/finite-element
    O16/O32/structured-bilinear
    O16/O32/structured-bicubic
    O16/O32/knn

    O16/H32/finite-element
    O16/H32/structured-bilinear
    O16/H32/structured-bicubic
    O16/H32/knn

    H16/O32/finite-element
    H16/O32/structured-bilinear
    H16/O32/structured-bicubic
    H16/O32/knn

    H16/H32/finite-element
    H16/H32/structured-bilinear
    H16/H32/structured-bicubic
    H16/H32/knn
  )

  unset(nprocs)
  list( APPEND nprocs 0 8 )

  foreach( case ${cases} )
    atlas_test_app_atlas_interpolations( ${case} "${nprocs}" )
  endforeach()

  # Due to currently unknown halo requirements for conservative interpolation, we only do this in serial here

  unset(cases)
  list( APPEND cases
    O16/O32/cons
    O16/O32/cons2

    O16/H32/cons
    O16/H32/cons2

    H16/O32/cons
    H16/O32/cons2

    H16/H32/cons
    H16/H32/cons2
  )

  unset(nprocs)
  list( APPEND nprocs 0 )

  foreach( case ${cases} )
    atlas_test_app_atlas_interpolations( ${case} "${nprocs}" )
  endforeach()



endif()