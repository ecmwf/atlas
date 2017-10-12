
if( HAVE_FCTEST )
    foreach(test
            fctest_griddistribution
            fctest_state)

        add_fctest( TARGET atlas_${test} SOURCES ${test}.F90 LINKER_LANGUAGE Fortran LIBS atlas_f )

    endforeach()
endif()

foreach(test
        test_domain
        test_field
        test_grid_ptr
        test_grids
        test_rotation
        test_state)

    ecbuild_add_test( TARGET atlas_${test} SOURCES ${test}.cc LIBS atlas )

endforeach()
