# (C) Copyright 1996-2017 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( MPI_CXX_FOUND AND MPI_Fortan_FOUND )

	ecbuild_info( "MPI" )

	ecbuild_info( "    MPI_CXX_COMPILER         : [${MPI_CXX_COMPILER}]" )
	ecbuild_info( "    MPI_CXX_INCLUDE_PATH     : [${MPI_CXX_INCLUDE_PATH}]" )
	ecbuild_info( "    MPI_CXX_LIBRARIES        : [${MPI_CXX_LIBRARIES}]" )

	ecbuild_info( "    MPI_Fortan_COMPILER      : [${MPI_Fortan_COMPILER}]" )
	ecbuild_info( "    MPI_Fortan_INCLUDE_PATH  : [${MPI_Fortan_INCLUDE_PATH}]" )
	ecbuild_info( "    MPI_Fortan_LIBRARIES     : [${MPI_Fortan_LIBRARIES}]" )

	ecbuild_info( "    MPIEXEC                  : [${MPIEXEC}]" )

endif()

if( CGAL_FOUND )

	ecbuild_info( "CGAL (${CGAL_VERSION})" )
	ecbuild_info( "    includes : [${CGAL_INCLUDE_DIRS}]" )
	ecbuild_info( "    libs     : [${CGAL_LIBRARY}]" )

endif()

if( ATLAS_HAVE_GRIDTOOLS_STORAGE )

    ecbuild_info( "GRIDTOOLS_STORAGE" )
    if( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST )
        ecbuild_info( "    BACKEND : [HOST]" )
    endif()    
    if( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA )
        ecbuild_info( "    BACKEND : [CUDA]" )
    endif()    

endif()

if( CUDA_FOUND )

    ecbuild_info( "CUDA (${CUDA_VERSION})" )
    ecbuild_info( "    CUDA_NVCC_COMPILER  : [${CUDA_NVCC_EXECUTABLE}]" )
    ecbuild_info( "    CUDA_CUDART_LIBRARY : [${CUDA_CUDART_LIBRARY}]" )
    ecbuild_info( "    CUDA_NVCC_FLAGS     : [${CUDA_NVCC_FLAGS}]" )

endif()

if( ATLAS_HAVE_ACC )

    ecbuild_info( "ACC" )
    ecbuild_info( "    ACC_C_COMPILER     : [${ACC_C_COMPILER}]" )
    ecbuild_info( "    ACC_C_FLAGS        : [${ACC_C_FLAGS}]" )
    ecbuild_info( "    ACC_Fortran_FLAGS  : [${ACC_Fortran_FLAGS}]" )

endif()