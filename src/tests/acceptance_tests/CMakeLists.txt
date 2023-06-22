# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable(
    TARGET atlas-atest-mgrids
    SOURCES atest_mgrids.cc
    LIBS atlas
    CONDITION atlas_HAVE_ATLAS_NUMERICS
)

if( HAVE_TESTS )
  ecbuild_add_option( FEATURE ACCEPTANCE_TESTS DEFAULT OFF )
endif()

if( HAVE_ACCEPTANCE_TESTS )

set( HAVE_ACCEPTANCE_TESTS_SMALL ON )
set( HAVE_ACCEPTANCE_TESTS_LARGE ON )

function( atlas_atest_mgrids category case nprocs )

    if( atlas_HAVE_ECTRANS )
      set( PARTITIONER "--partitioner=trans" )
    endif()

    string(REPLACE "/" ";" case ${case} )
    list( GET case 0 gridA )
    list( GET case 1 gridB )
    unset( checksums )
    unset( case_tests )
    foreach( n ${nprocs} )
        if( MPI_SLOTS GREATER_EQUAL ${n} )
            foreach( t 1 4 12 )
                set( test_name  atlas_atest_${category}_mgrids_${gridA}_to_${gridB}_n${n}_t${t} )
                ecbuild_add_test( TARGET ${test_name}
                    COMMAND atlas-atest-mgrids
                    ARGS    --gridA=${gridA}
                            --gridB=${gridB}
                            --haloA=3
                            --haloB=3
                            --display-name=${test_name}
                            --checksum=${test_name}.checksum
                            --matrix-free
                            ${PARTITIONER}
                    MPI ${n}
                    OMP ${t}
                    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                )
                list( APPEND checksums ${test_name}.checksum )
                list( APPEND case_tests ${test_name} )
            endforeach()
        endif()
    endforeach()
    add_test(
        NAME atlas_atest_${category}_mgrids_${gridA}_to_${gridB}_BitReproducible
        COMMAND ${CMAKE_COMMAND}
        "-DFILES=${checksums}" -P ${CMAKE_CURRENT_SOURCE_DIR}/compare.cmake
    )
    foreach( dep ${case_tests} )
        set_tests_properties( atlas_atest_${category}_mgrids_${gridA}_to_${gridB}_BitReproducible PROPERTIES DEPENDS ${dep} )
    endforeach()
endfunction()

if( atlas_HAVE_ATLAS_NUMERICS )
if( HAVE_ACCEPTANCE_TESTS_SMALL )
    unset( cases )
    list( APPEND cases
        O32/O16
        O32/L16
        O64/F20
    )
    list( APPEND nprocs
        1
        4
        32
    )
    foreach( case ${cases} )
        atlas_atest_mgrids( small ${case} "${nprocs}" )
    endforeach()
endif()

if( HAVE_ACCEPTANCE_TESTS_LARGE )
    unset( cases )
    list( APPEND cases
        O640/O320
    )
    set( nprocs 1 4 36 72 )
    foreach( case ${cases} )
        atlas_atest_mgrids( large ${case} "${nprocs}" )
    endforeach()
endif()
endif()
endif()

