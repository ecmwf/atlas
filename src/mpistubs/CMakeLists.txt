# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_library( TARGET mpistubs
  AUTO_VERSION
  INSTALL_HEADERS ALL
  HEADER_DESTINATION  ${INSTALL_INCLUDE_DIR}/mpistubs
  SOURCES   mpistubs.h
            mpistubs.c )

ecbuild_add_library( TARGET mpistubs_f
  AUTO_VERSION
  CONDITION HAVE_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED
  SOURCES   mpistubs.F90 )
