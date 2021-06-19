set( ENABLE_TRANS               ON  CACHE BOOL "Enable TRANS" )
set( ENABLE_BOUNDSCHECKING      ON  CACHE BOOL "Enable bounds checking")
set( ENABLE_TESSELATION         OFF CACHE BOOL "Disable CGAL" ) # cgal is old in leap42
set( ENABLE_OMP_CXX             OFF CACHE BOOL "Disable OpenMP for C++" ) # because of problems with clang OpenMP
