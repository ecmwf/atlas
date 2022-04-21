### OMP ...

if( ENABLE_OMP OR NOT DEFINED ENABLE_OMP )
    find_package( OpenMP COMPONENTS CXX ${Fortran} )
endif()

ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "support for OpenMP shared memory parallelism"
                    CONDITION OpenMP_Fortran_FOUND OR OpenMP_CXX_FOUND )
ecbuild_add_option( FEATURE OMP_Fortran
                    DESCRIPTION "support for Fortran OpenMP shared memory parallelism"
                    CONDITION HAVE_OMP AND OpenMP_Fortran_FOUND )

ecbuild_add_option( FEATURE OMP_CXX
                    DESCRIPTION "support for CXX OpenMP shared memory parallelism"
                    CONDITION HAVE_OMP AND OpenMP_CXX_FOUND )
if( TARGET OpenMP::OpenMP_CXX )
  set( OMP_CXX OpenMP::OpenMP_CXX )
endif()
if( TARGET OpenMP::OpenMP_Fortran )
  set( OMP_Fortran OpenMP::OpenMP_Fortran )
endif()

if( HAVE_OMP_CXX )

  if( NOT CMAKE_CXX_COMPILER_ID MATCHES Clang )
    set( ATLAS_OMP_TASK_SUPPORTED 1 )
  endif()

  if( NOT DEFINED ATLAS_OMP_TASK_SUPPORTED )
    try_run( execute_result compile_result
                     ${CMAKE_CURRENT_BINARY_DIR}
                     ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_task.cc
                     LINK_LIBRARIES ${OMP_CXX}
                     COMPILE_OUTPUT_VARIABLE compile_output
                     RUN_OUTPUT_VARIABLE execute_output )

    ecbuild_debug("Compiling and running ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_task.cc")
    ecbuild_debug_var( compile_result )
    ecbuild_debug_var( compile_output )
    ecbuild_debug_var( execute_result )
    ecbuild_debug_var( execute_output )

    if( compile_result )
      if( execute_result MATCHES 0 )
        set( ATLAS_OMP_TASK_SUPPORTED 1 )
      else()
        ecbuild_info("    Compiler failed to correctly run program with 'omp task' pragma."
                     "Sorting with OMP is disabled.")
        set( ATLAS_OMP_TASK_SUPPORTED 0 )
      endif()
    else()
      set( ATLAS_OMP_TASK_SUPPORTED 0 )
    endif()
  endif()


  if( ATLAS_OMP_TASK_SUPPORTED )
    if( NOT CMAKE_CXX_COMPILER_ID MATCHES Clang )
      set( ATLAS_OMP_TASK_UNTIED_SUPPORTED 1 )
    endif()

    if( NOT DEFINED ATLAS_OMP_TASK_UNTIED_SUPPORTED )
      try_run( execute_result compile_result
                       ${CMAKE_CURRENT_BINARY_DIR}
                       ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_untied.cc
                       LINK_LIBRARIES ${OMP_CXX}
                       COMPILE_OUTPUT_VARIABLE compile_output
                       RUN_OUTPUT_VARIABLE execute_output )

      ecbuild_debug("Compiling and running ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_untied.cc")
      ecbuild_debug_var( compile_result )
      ecbuild_debug_var( compile_output )
      ecbuild_debug_var( execute_result )
      ecbuild_debug_var( execute_output )

      if( compile_result )
        if( execute_result MATCHES 0 )
          set( ATLAS_OMP_TASK_UNTIED_SUPPORTED 1 )
        else()
          ecbuild_info("    Compiler failed to run program with omp pragma with 'untied if' construct."
                       "Workaround will be enabled.")
          set( ATLAS_OMP_TASK_UNTIED_SUPPORTED 0 )
        endif()
      else()
        set( ATLAS_OMP_TASK_UNTIED_SUPPORTED 0 )
      endif()
    endif()
  else()
    set( ATLAS_OMP_TASK_UNTIED_SUPPORTED 0 )
  endif()
else()
  set( ATLAS_OMP_TASK_SUPPORTED 0 )
  set( ATLAS_OMP_TASK_UNTIED_SUPPORTED 0 )
endif()
