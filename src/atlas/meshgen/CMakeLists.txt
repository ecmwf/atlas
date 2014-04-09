list( APPEND atlas_meshgen_src
  RGG.cpp
  RGG.hpp
  T159.cpp
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
