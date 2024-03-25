
ecbuild_add_library( TARGET atlas_io

    INSTALL_HEADERS    ALL
    HEADER_DESTINATION include/atlas_io
    PUBLIC_LIBS        eckit_codec
    PUBLIC_INCLUDES
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/eckit_codec_adaptor/src>

    SOURCES
        atlas-io.h
        Trace.cc
        Trace.h
        detail/BlackMagic.h
    LINKER_LANGUAGE CXX
)

target_compile_features( atlas_io PUBLIC cxx_std_17 )

