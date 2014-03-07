ecbuild_add_test(
        TARGET      cgal_mesh_gen_on_surface
        SOURCES     cgal_mesh_gen_on_surface.cc
        CONDITION   HAVE_TRIANGULATION AND CGAL_FOUND
	INCLUDES    ${CGAL_INCLUDE_DIRS}
        LIBS        atlas ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES}
	)

ecbuild_add_test(
        TARGET      cgal_mesh_gen_from_points
        SOURCES     cgal_mesh_gen_from_points.cc
        CONDITION   HAVE_TRIANGULATION AND CGAL_FOUND
        INCLUDES    ${CGAL_INCLUDE_DIRS}
        LIBS        atlas ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES}
        )
