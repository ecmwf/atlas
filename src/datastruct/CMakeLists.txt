
function(generate_fortran_bindings output filename)
    get_filename_component(base ${filename} NAME_WE)
    set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
    set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}.f90) 
    set(${output} ${${output}} ${outfile} PARENT_SCOPE)
    add_custom_command(
        OUTPUT ${outfile}
        COMMAND python ${CMAKE_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
        DEPENDS ${filename} )
    set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE)
endfunction()

generate_fortran_bindings(FORTRAN_BINDINGS Field.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS FieldSet.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS FunctionSpace.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS Mesh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS Metadata.hpp)


add_library( datastruct
  # C++ part
    Field.hpp 
    Field.cpp
    FieldSet.hpp 
    FieldSet.cpp 
    FunctionSpace.hpp 
    FunctionSpace.cpp
    Mesh.hpp
    Mesh.cpp
    Metadata.hpp
    Metadata.cpp
    Comm.hpp
    Comm.cpp

  # Fortran part
    ${FORTRAN_BINDINGS}
    datastruct.F90
) 
target_link_libraries( datastruct common ${TPL} )

