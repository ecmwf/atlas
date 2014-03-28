ecbuild_add_test(
        TARGET      cgal_mesh_gen_from_points
        SOURCES     cgal_mesh_gen_from_points.cc
		CONDITION   HAVE_TESSELATION AND CGAL_FOUND
		LIBS        atlas_grid atlas )

