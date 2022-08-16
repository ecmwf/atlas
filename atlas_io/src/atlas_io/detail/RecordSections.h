/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstdint>
#include <string>

#include "eckit/types/FixedString.h"

#include "atlas_io/detail/Endian.h"
#include "atlas_io/detail/Time.h"
#include "atlas_io/detail/Version.h"

namespace atlas {
namespace io {

// ------------------------------------------------------------------------------------------------------------------------------------

struct RecordBegin {
    eckit::FixedString<8> string{"ATLAS-IO"};  ///< 16 ATLAS-IO string
    eckit::FixedString<8> spare{"\n"};         ///<

    bool valid() const { return string == "ATLAS-IO"; }
    std::string str() const { return string; }
};

struct RecordEnd {  // 32 bytes
    static constexpr size_t bytes = 32;

    eckit::FixedString<1> eol{"\n"};
    eckit::FixedString<12> string{"ATLAS-IO-END"};
    eckit::FixedString<19> padding{"               \n\n\n\n"};

    bool valid() const { return string == "ATLAS-IO-END"; }
    std::string str() const { return string; }
};

// ------------------------------------------------------------------------------------------------------------------------------------

struct RecordHead {
    static constexpr size_t bytes    = 256;
    static constexpr size_t padding_ = bytes - 16 - 8 - 16 - 8 * 4 - 64 - 8 * 2 - 4 - 1;

    RecordBegin begin;                              ///< 16 beginning of record
    Version version;                                ///< 8  version of this record
    Time time;                                      ///< 16 time since system_clock epoch (1970-1-1 00:00)
    uint64_t record_length{0};                      ///<  8 length of entire record
    eckit::FixedString<8> metadata_format{"yaml"};  ///<  8 format of metadata section in this record
    std::uint64_t metadata_offset{bytes};           ///<  8 offset where metadata section starts
    std::uint64_t metadata_length{0};               ///<  8 length of metadata section
    eckit::FixedString<64> metadata_checksum;       ///< 64 checksum of metadata
    std::uint64_t index_offset{0};                  ///<  8  offset where data section starts
    std::uint64_t index_length;                     ///<  8 length of data section
    std::uint32_t magic_number{1234};        ///<  4 number 1234 encoded in binary, used to detect encoded endianness
    eckit::FixedString<padding_> padding__;  ///<  Extra padding to get to <bytes>
    eckit::FixedString<1> eol{"\n"};

    static constexpr size_t size() { return bytes; }  ///< Size in bytes of this section
    Endian endian() const;                            ///< Endianness determined from magic_number
    bool valid() const { return begin.valid(); }      ///< Check if this is a valid RecordRoot
    std::string str() const { return begin.str(); }
};

// ------------------------------------------------------------------------------------------------------------------------------------

struct RecordMetadataSection {
    struct Begin {  // 32 bytes
        static constexpr size_t bytes = 32;

        eckit::FixedString<1> eol{"\n"};
        eckit::FixedString<14> string{"METADATA-BEGIN"};
        eckit::FixedString<17> padding{"                \n"};

        bool valid() const { return string == "METADATA-BEGIN"; }
        std::string str() const { return string; }
    };

    struct End {  // 32 bytes
        static constexpr size_t bytes = 32;

        eckit::FixedString<1> eol{"\n"};
        eckit::FixedString<12> string{"METADATA-END"};
        eckit::FixedString<19> padding{"                  \n"};

        bool valid() const { return string == "METADATA-END"; }
        std::string str() const { return string; }
    };
};

// ------------------------------------------------------------------------------------------------------------------------------------

struct RecordDataIndexSection {
    struct Begin {  // 32 bytes
        static constexpr size_t bytes = 32;

        eckit::FixedString<1> eol{"\n"};
        eckit::FixedString<11> string{"INDEX-BEGIN"};
        eckit::FixedString<20> padding{"                   \n"};

        bool valid() const { return string == "INDEX-BEGIN"; }
        std::string str() const { return string; }
    };

    struct End {  // 32 bytes
        static constexpr size_t bytes = 32;

        eckit::FixedString<1> eol{"\n"};
        eckit::FixedString<9> string{"INDEX-END"};
        eckit::FixedString<22> padding{"                     \n"};

        bool valid() const { return string == "INDEX-END"; }
        std::string str() const { return string; }
    };

    struct Entry {
        static constexpr size_t bytes = 80;

        std::uint64_t offset;
        std::uint64_t length;
        eckit::FixedString<64> checksum;
    };
};

// ------------------------------------------------------------------------------------------------------------------------------------

struct RecordDataSection {
    struct Begin {  // 32 bytes
        static constexpr size_t bytes = 32;

        eckit::FixedString<1> eol{"\n"};
        eckit::FixedString<10> string{"DATA-BEGIN"};
        eckit::FixedString<21> padding{"                    \n"};

        bool valid() const { return string == "DATA-BEGIN"; }
        std::string str() const { return string; }
    };
    struct End {  // 32 bytes
        static constexpr size_t bytes = 32;

        eckit::FixedString<1> eol{"\n"};
        eckit::FixedString<8> string{"DATA-END"};
        eckit::FixedString<23> padding{"                      \n"};

        bool valid() const { return string == "DATA-END"; }
        std::string str() const { return string; }
    };
};

// ------------------------------------------------------------------------------------------------------------------------------------


}  // namespace io
}  // namespace atlas
