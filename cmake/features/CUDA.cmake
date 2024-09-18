
# ecbuild_add_option( FEATURE CUDA
#                     DESCRIPTION "Enable CUDA support"
#                     DEFAULT OFF
#                   )
# ecbuild_add_option( FEATURE HIP
#                     DESCRIPTION "Enable CUDA support"
#                     DEFAULT OFF
#                   )

set( atlas_HAVE_CUDA 0 )
set( atlas_HAVE_HIP  0 )
set( atlas_HAVE_GPU  0 )

if( hic_HAVE_CUDA )
  enable_language( CUDA )
  ecbuild_info( "CUDA language enabled" )
  find_package(CUDAToolkit REQUIRED)
  set( atlas_HAVE_CUDA 1 )
  set( atlas_HAVE_GPU  1 )
elseif( hic_HAVE_HIP )
  enable_language( HIP )
  ecbuild_info( "HIP language enabled" )
  find_package(hip CONFIG REQUIRED)
  set( atlas_HAVE_HIP 1 )
  set( atlas_HAVE_GPU 1 )
endif()

set( HAVE_CUDA ${atlas_HAVE_CUDA} )
set( HAVE_HIP  ${atlas_HAVE_HIP} )
set( HAVE_GPU  ${atlas_HAVE_GPU} )

if( HAVE_GPU )
  ecbuild_info("GPU support enabled")
else()
  ecbuild_info("GPU support not enabled")
endif()