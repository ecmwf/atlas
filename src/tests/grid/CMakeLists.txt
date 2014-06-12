ecbuild_add_test( TARGET      test_RegularLatLonGrid
			   BOOST
               SOURCES     test_RegularLatLonGrid.cc
               CONDITION   HAVE_GRID AND HAVE_GRIB
               LIBS        atlas_grid  
                )

ecbuild_add_test(	TARGET      test_latlon
					BOOST
					SOURCES     test_latlon.cc
					CONDITION   HAVE_GRID
					LIBS        atlas_grid  )
            
ecbuild_add_test(   TARGET      test_gaussian
					BOOST
					SOURCES     test_gaussian.cc
					CONDITION	HAVE_GRID
					LIBS        atlas_grid )

ecbuild_add_test(   TARGET      test_field
					BOOST
					SOURCES     test_field.cc
          CONDITION	HAVE_GRID AND HAVE_GRIB
					LIBS        atlas_grid )
