
### MPI ...
if( NOT eckit_HAVE_MPI )
  ecbuild_warn("ecKit has been compiled without MPI. This causes Atlas to not be able to run parallel jobs.")
  set( atlas_HAVE_MPI 0 )
else()
  set( atlas_HAVE_MPI 1 )
endif()
