if( NOT HAVE_MPI )

	include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/mpi )
	list( APPEND TPL mpi )

	add_subdirectory( mpi)

endif()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_BINARY_DIR} )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/atlas )

add_subdirectory( atlas )
add_subdirectory( tests )
