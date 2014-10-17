ecbuild_add_test( TARGET      atlas_test_Grid
			      BOOST
                  SOURCES     test_Grid.cc
				  CONDITION   ECKIT_HAVE_GRIB AND OFF
				  LIBS        atlas
                )

ecbuild_add_test( TARGET      atlas_test_GridSpec
                  BOOST
                  SOURCES     test_GridSpec.cc
				  CONDITION   ECKIT_HAVE_GRIB AND OFF
				  LIBS        atlas
                )


#ecbuild_add_test(	TARGET      test_latlon
#					BOOST
#					SOURCES     test_latlon.cc
#					CONDITION   HAVE_GRID
#					LIBS        atlas  )

#ecbuild_add_test(   TARGET      test_gaussian
#					BOOST
#					SOURCES     test_gaussian.cc
#					CONDITION	HAVE_GRID
#					LIBS        atlas )

ecbuild_add_test(   TARGET      atlas_test_field
					BOOST
					SOURCES     test_field.cc
					CONDITION	ECKIT_HAVE_GRIB
					LIBS        atlas )

list( APPEND test_grid_from_grib_data
CMC_hrdps_east_TMP_TGL_2_ps2.5km_2014081512_P000-00.grib2
ll005005.grib
ll0101.grib
ll11.grib
ll55.grib
n1024.grib
#n1280.grib
n160.grib
n256.grib
n320.grib
n40.grib
n400.grib
n512.grib
n60.grib
rgg_n640.grib )

foreach( f ${test_grid_from_grib_data} )
  set(tgt data_${f})
  ecbuild_get_test_data( TARGET "${tgt}" NAME "${f}" DIRNAME ${PROJECT_NAME}/grid NOCHECK )
  list(APPEND tgts "${tgt}")
endforeach()

add_custom_target( atlas_tests_grid_get_gribs DEPENDS ${tgts} )

ecbuild_add_test(   TARGET      atlas_test_grid_from_grib
					BOOST
					SOURCES     test_grid_from_grib.cc
                    DEPENDS     atlas_tests_grid_get_gribs
					CONDITION	ECKIT_HAVE_GRIB AND OFF
					LIBS        atlas )
