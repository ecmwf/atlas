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
					LIBS        atlas
                    MPI         1 )
