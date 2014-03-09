ecbuild_add_test(
        TARGET      cgal_mesh_gen_on_surface
        SOURCES     cgal_mesh_gen_on_surface.cc
		CONDITION   HAVE_TRIANGULATION AND CGAL_FOUND AND OFF
		INCLUDES    ${CGAL_INCLUDE_DIRS}
        LIBS        atlas ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES}
	)

ecbuild_add_test(
        TARGET      cgal_mesh_gen_from_points
        SOURCES     cgal_mesh_gen_from_points.cc
		CONDITION   HAVE_TRIANGULATION AND CGAL_FOUND
		LIBS        atlas )

ecbuild_add_test(
		TARGET      cgal_mesh_gen_from_grib
		SOURCES     cgal_mesh_gen_from_grib.cc
		CONDITION   HAVE_TRIANGULATION AND CGAL_FOUND AND GRIB_API_FOUND
		INCLUDES    ${GRIB_API_INCLUDE_DIRS}
		LIBS        atlas ${GRIB_API_LIBRARIES} )
