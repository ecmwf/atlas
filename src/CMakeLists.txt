### mpistubs

if( NOT HAVE_MPI )
  include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/mpistubs )
  list( APPEND TPL mpistubs )
  add_subdirectory( mpistubs )
endif()

### sources

add_subdirectory( atlas )

add_subdirectory( apps )

add_subdirectory( tests )

add_subdirectory( sandbox )
