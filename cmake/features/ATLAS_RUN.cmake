### use of atlas-run for tests
ecbuild_add_option( FEATURE ATLAS_RUN
                    DEFAULT ON
                    DESCRIPTION "Use atlas/tools/atlas-run to run atlas tests" )

if( HAVE_ATLAS_RUN )
  set( MPIEXEC_EXECUTABLE ${CMAKE_CURRENT_SOURCE_DIR}/tools/atlas-run )
  set( MPIEXEC_NUMPROC_FLAG='-n' )
  set( CMAKE_CROSSCOMPILING_EMULATOR ${CMAKE_CURRENT_SOURCE_DIR}/tools/atlas-run )
endif()

