
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
generate_fortran_bindings(FORTRAN_BINDINGS HaloExchange.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS Gather.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS Gmsh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS BuildPeriodicBoundaries.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS BuildEdges.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS BuildDualMesh.hpp)


add_library( atlas
  # C++ part
    BuildDualMesh.hpp
    BuildDualMesh.cpp
    BuildEdges.hpp
    BuildEdges.cpp
    BuildHalo.hpp
    BuildHalo.cpp
    BuildPeriodicBoundaries.hpp
    BuildPeriodicBoundaries.cpp
    Field.hpp 
    Field.cpp
    FieldSet.hpp 
    FieldSet.cpp 
    FunctionSpace.hpp 
    FunctionSpace.cpp
    Gmsh.hpp
    Gmsh.cpp
    Mesh.hpp
    Mesh.cpp
    Metadata.hpp
    Metadata.cpp
    Parameters.hpp
    Partitioner.hpp
    Partitioner.cpp
    HaloExchange.hpp
    HaloExchange.cpp
    Gather.hpp
    Gather.cpp
    MPL.hpp

  # Fortran part
    ${FORTRAN_BINDINGS}
    atlas.F90
) 
target_link_libraries( atlas ${TPL} )

