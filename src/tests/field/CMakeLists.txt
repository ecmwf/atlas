# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_field
    CONDITION       ATLAS_HAVE_GRIDTOOLS_STORAGE AND ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_field.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_field_gpu
    CONDITION       ATLAS_HAVE_GRIDTOOLS_STORAGE AND ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_field_gpu.F90 external_acc_routine.F90
    LIBS            atlas_f
  )
  if( TARGET atlas_fctest_field_gpu )
    target_compile_options( atlas_fctest_field_gpu PUBLIC ${ACC_Fortran_FLAGS} )
    target_link_libraries( atlas_fctest_field_gpu ${ACC_Fortran_FLAGS} )
  endif()

endif()

