# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( atlas_HAVE_ACC AND OpenACC_Fortran_FOUND )
    ecbuild_add_executable(
        TARGET  atlas-acc-fields
        SOURCES atlas-acc-fields.F90
        LIBS    atlas_f OpenACC::OpenACC_Fortran
        LINKER_LANGUAGE Fortran
        NOINSTALL
    )
endif()
