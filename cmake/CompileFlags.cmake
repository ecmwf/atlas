if( CMAKE_CXX_COMPILER_ID MATCHES Cray )

  ecbuild_add_cxx_flags("-hnomessage=3140") # colon separated numbers
  ecbuild_add_fortran_flags("-hnomessage=3140") # colon separated numbers

# CC-3140 crayc++: WARNING File = atlas/functionspace/NodeColumns.cc, Line = 1, Column = 1
#  The IPA optimization level was changed to "1" due to the presence of OMP
#          directives, ACC directives, or ASM intrinsics.

endif()
