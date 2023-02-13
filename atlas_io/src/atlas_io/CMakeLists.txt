configure_file( detail/defines.h.in   detail/defines.h   )
install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/detail/defines.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/atlas_io/detail
)

ecbuild_add_library( TARGET atlas_io

    INSTALL_HEADERS    ALL
    HEADER_DESTINATION include/atlas_io
    PUBLIC_LIBS        eckit
    PUBLIC_INCLUDES
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>

    SOURCES
        atlas-io.h
        Data.cc
        Data.h
        detail/Assert.h
        detail/Base64.cc
        detail/Base64.h
        detail/Checksum.h
        detail/Checksum.cc
        detail/DataInfo.h
        detail/DataType.cc
        detail/DataType.h
        detail/Decoder.cc
        detail/Decoder.h
        detail/Defaults.h
        detail/Encoder.cc
        detail/Encoder.h
        detail/Endian.h
        detail/Link.cc
        detail/Link.h
        detail/ParsedRecord.h
        detail/RecordInfo.h
        detail/RecordSections.h
        detail/Reference.h
        detail/sfinae.h
        detail/StaticAssert.h
        detail/tag.h
        detail/Time.cc
        detail/Time.h
        detail/Type.h
        detail/TypeTraits.h
        detail/Version.h
        Exceptions.cc
        Exceptions.h
        FileStream.cc
        FileStream.h
        Metadata.cc
        Metadata.h
        print/TableFormat.cc
        print/TableFormat.h
        print/JSONFormat.cc
        print/JSONFormat.h
        print/Bytes.cc
        print/Bytes.h
        ReadRequest.cc
        ReadRequest.h
        Record.cc
        Record.h
        RecordItem.cc
        RecordItem.h
        RecordItemReader.cc
        RecordItemReader.h
        RecordPrinter.cc
        RecordPrinter.h
        RecordReader.cc
        RecordReader.h
        RecordWriter.cc
        RecordWriter.h
        Session.cc
        Session.h
        Stream.cc
        Stream.h
        Trace.cc
        Trace.h
        types/array.h
        types/array/ArrayMetadata.cc
        types/array/ArrayMetadata.h
        types/array/ArrayReference.cc
        types/array/ArrayReference.h
        types/array/adaptors/StdArrayAdaptor.h
        types/array/adaptors/StdVectorAdaptor.h
        types/string.h
        types/scalar.h
        types/scalar.cc
        ${CMAKE_CURRENT_BINARY_DIR}/detail/defines.h
)

target_compile_features( atlas_io PUBLIC cxx_std_17 )
