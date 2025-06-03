function( add_sandbox target )
  set( options           "")
  set( single_value_args "")
  set( multi_value_args  "")
  cmake_parse_arguments(PARSE_ARGV 1 _PAR "${options}" "${single_value_args}" "${multi_value_args}")
  set( sources ${_PAR_UNPARSED_ARGUMENTS})
  ecbuild_add_executable(TARGET ${target} SOURCES ${sources} LIBS pluto NOINSTALL)
  pluto_declare_hic_source_files( ${sources} )
  pluto_target_linker_language_hic( ${target} )
endfunction()


add_sandbox(pluto-sandbox sandbox.cc kernel.cc)
