# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

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
