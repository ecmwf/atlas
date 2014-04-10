list( APPEND atlas_meshgen_src
  EqualAreaPartitioner.cpp
  EqualAreaPartitioner.hpp
  RGG.cpp
  RGG.hpp
  T63.cpp
  T95.cpp
  T159.cpp
  T255.cpp
  T511.cpp
  T1279.cpp
  T2047.cpp
)

ecbuild_add_library(TARGET atlas_meshgen
					AUTO_VERSION
					INSTALL_HEADERS ALL
					SOURCES
						${atlas_meshgen_src}
					HEADER_DESTINATION
						${INSTALL_INCLUDE_DIR}/atlas/meshgen
					LIBS
						atlas )
