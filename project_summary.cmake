if( MPI_CXX_FOUND AND MPI_Fortan_FOUND )

	message( STATUS "MPI")

	message( STATUS "	MPI_CXX_COMPILER      : [${MPI_CXX_COMPILER}]")
	message( STATUS "	MPI_CXX_INCLUDE_PATH  : [${MPI_CXX_INCLUDE_PATH}]")
	message( STATUS "	MPI_CXX_LIBRARIES     : [${MPI_CXX_LIBRARIES}]")

	message( STATUS "	MPI_Fortan_COMPILER      : [${MPI_Fortan_COMPILER}]")
	message( STATUS "	MPI_Fortan_INCLUDE_PATH  : [${MPI_Fortan_INCLUDE_PATH}]")
	message( STATUS "	MPI_Fortan_LIBRARIES     : [${MPI_Fortan_LIBRARIES}]")

	message( STATUS "	MPIEXEC               : [${MPIEXEC}]")

endif()

if( CGAL_FOUND )

	message( STATUS "CGAL (${CGAL_VERSION})")
	message( STATUS "	includes : [${CGAL_INCLUDE_DIRS}]")
	message( STATUS "	libs     : [${CGAL_LIBRARY}]")

endif()