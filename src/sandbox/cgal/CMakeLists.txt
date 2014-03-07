ecbuild_add_test(
	TARGET      cgal_mesh_gen
	SOURCES     cgal_mesh_gen.cc
	CONDITION   HAVE_TRIANGULATION AND CGAL_FOUND
	INCLUDES    ${CGAL_INCLUDE_DIRS}
	LIBS        atlas ${CGAL_LIBRARY}
	)
