list( APPEND atlas_mpi_srcs
        Comm.h
        Comm.cc
)

ecbuild_add_library( TARGET atlas_mpi
                                         INSTALL_HEADERS ALL
                                         HEADER_DESTINATION
                                                ${INSTALL_INCLUDE_DIR}/atlas/mpi
                                         CONDITION HAVE_MPI AND HAVE_SANDBOX
                                         SOURCES
                                                ${atlas_mpi_srcs}
                                         INCLUDES
                                                ${MPI_CXX_INCLUDE_PATH}
                                         LIBS
                                                atlas ${MPI_CXX_LIBRARIES} )
ecbuild_add_test(
	TARGET      mpi_comm
	SOURCES     mpi_comm.cc 
	CONDITION   HAVE_MPI
	INCLUDES    ${MPI_CXX_INCLUDE_PATH} 
	LIBS        atlas_mpi atlas ${MPI_CXX_LIBRARIES} 
	)

