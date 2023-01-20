/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas_io/Record.h"

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/filesystem/URI.h"

#include "atlas_io/Exceptions.h"
#include "atlas_io/Trace.h"
#include "atlas_io/detail/Assert.h"
#include "atlas_io/detail/ParsedRecord.h"
#include "atlas_io/detail/Version.h"

namespace atlas {
namespace io {

namespace {

//---------------------------------------------------------------------------------------------------------------------

template <typename IStream, typename Struct>
inline size_t read_struct(IStream& in, Struct& s) {
    static_assert(Struct::bytes == sizeof(Struct), "");
    return in.read(reinterpret_cast<char*>(&s), sizeof(Struct));
}

//---------------------------------------------------------------------------------------------------------------------

template <typename Struct, typename IStream>
inline Struct read_struct(IStream& in) {
    Struct s;
    if (read_struct(in, s) != sizeof(Struct)) {
        throw InvalidRecord("Unexpected EOF reached");
    }
    return s;
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace

//---------------------------------------------------------------------------------------------------------------------

Endian RecordHead::endian() const {
    if (magic_number == 1234) {
        return Endian::native;
    }
    else if (magic_number == 3523477504) {
        return Endian::swapped;
    }
    throw Exception("Mixed endianness is not supported", Here());
}

//---------------------------------------------------------------------------------------------------------------------

std::string Record::URI::str() const {
    eckit::URI uri("file", eckit::PathName(path));
    uri.query("offset", std::to_string(offset));
    return uri.asRawString();
}

//---------------------------------------------------------------------------------------------------------------------

Record::Record(): record_(new ParsedRecord()) {}

Record::Record(const Record& other) = default;

//---------------------------------------------------------------------------------------------------------------------

bool Record::empty() const {
    return record_->head.record_length == 0;
}

//---------------------------------------------------------------------------------------------------------------------

const Metadata& Record::metadata(const std::string& key) const {
    if (record_->items.find(key) == record_->items.end()) {
        throw Exception("Record does not contain key \"" + key + "\"", Here());
    }
    return record_->items.at(key);
}

//---------------------------------------------------------------------------------------------------------------------

Endian Record::endian() const {
    return record_->head.endian();
}

//---------------------------------------------------------------------------------------------------------------------

Version Record::version() const {
    return record_->head.version;
}

//---------------------------------------------------------------------------------------------------------------------

Time Record::time() const {
    return record_->head.time;
}

//---------------------------------------------------------------------------------------------------------------------

uint64_t Record::size() const {
    return record_->head.record_length;
}

//---------------------------------------------------------------------------------------------------------------------

const std::vector<std::string>& Record::keys() const {
    return record_->keys;
}

//---------------------------------------------------------------------------------------------------------------------

bool Record::has(const std::string& key) {
    return record_->items.find(key) != record_->items.end();
}

//---------------------------------------------------------------------------------------------------------------------

Record::operator const ParsedRecord&() const {
    return *record_;
}

//---------------------------------------------------------------------------------------------------------------------

static void parse_record(ParsedRecord& record, const std::string& key, const Metadata& metadata) {
    if (metadata.type() || metadata.link()) {
        record.items.emplace(key, metadata);
        record.keys.emplace_back(key);
    }
    else {
        for (auto& next_key : metadata.keys()) {
            parse_record(record, key + "." + next_key, metadata.getSubConfiguration(next_key));
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

Record& Record::read(Stream& in, bool read_to_end) {
    if (not empty()) {
        return *this;
    }

    ATLAS_IO_TRACE("read_metadata");
    auto& r = record_->head;

    auto rbegin = in.position();

    // Begin Record
    // ------------
    if (atlas::io::read_struct(in, r) != sizeof(r)) {
        if (in.position() > sizeof(r.begin.string)) {
            if (not r.begin.valid()) {
                std::stringstream err;
                err << "Format is not recognized. Received: " << r.begin.string;
                throw InvalidRecord(err.str());
            }
        }
        throw InvalidRecord("Unexpected EOF reached");
    }
    if (not r.valid()) {
        std::stringstream err;
        err << "Format is not recognized. Received: " << r.begin.string;
        throw InvalidRecord(err.str());
    }
    if (r.version < RecordHead{}.version) {
        throw InvalidRecord("Version of record (" + r.version.str() + ") is too old.");
    }

    if (r.metadata_length < sizeof(RecordMetadataSection::Begin) + sizeof(RecordMetadataSection::End)) {
        throw InvalidRecord("Unexpected metadata section length: " + std::to_string(r.metadata_length) + " < " +
                            std::to_string(sizeof(RecordMetadataSection::Begin) + sizeof(RecordMetadataSection::End)));
    }

    if (r.index_length < sizeof(RecordDataIndexSection::Begin) + sizeof(RecordDataIndexSection::End)) {
        throw InvalidRecord("Unexpected data index section length.");
    }

    r.metadata_offset += rbegin;
    r.index_offset += rbegin;

    // Metadata section
    // ----------------
    in.seek(r.metadata_offset);
    auto metadata_begin = atlas::io::read_struct<RecordMetadataSection::Begin>(in);
    if (not metadata_begin.valid()) {
        throw InvalidRecord("Metadata section is not valid. Invalid section begin marker: [" + metadata_begin.str() +
                            "]");
    }
    std::string metadata_str;
    metadata_str.resize(size_t(r.metadata_length) - sizeof(RecordMetadataSection::Begin) -
                        sizeof(RecordMetadataSection::End));
    if (in.read(const_cast<char*>(metadata_str.data()), metadata_str.size()) != metadata_str.size()) {
        throw InvalidRecord("Unexpected EOF reached");
    }
    auto metadata_end = atlas::io::read_struct<RecordMetadataSection::End>(in);
    if (not metadata_end.valid()) {
        throw InvalidRecord("Metadata section is not valid. Invalid section end marker: [" + metadata_end.str() + "]");
    }
    Checksum encoded_metadata_checksum(r.metadata_checksum);
    Checksum computed_metadata_checksum(
        atlas::io::checksum(metadata_str.data(), metadata_str.size(), encoded_metadata_checksum.algorithm()));
    if (computed_metadata_checksum.available() && encoded_metadata_checksum.str() != computed_metadata_checksum.str()) {
        std::stringstream err;
        err << "Mismatch in metadata checksum.\n";
        err << "        Encoded:  [" << encoded_metadata_checksum.str() << "].\n";
        err << "        Computed: [" << computed_metadata_checksum.str() << "].";
        throw DataCorruption(err.str());
    }

    ATLAS_IO_ASSERT(r.metadata_format == "yaml");
    Metadata metadata = eckit::YAMLConfiguration(metadata_str);

    for (auto& key : metadata.keys()) {
        parse_record(*record_, key, metadata.getSubConfiguration(key));
    }


    // DataIndex section
    // -----------------
    in.seek(r.index_offset);
    auto index_begin = atlas::io::read_struct<RecordDataIndexSection::Begin>(in);
    if (not index_begin.valid()) {
        throw InvalidRecord("Data index section is not valid. Invalid section begin marker: [" + index_begin.str() +
                            "]");
    }
    const auto index_length =
        (size_t(r.index_length) - sizeof(RecordDataIndexSection::Begin) - sizeof(RecordDataIndexSection::End));
    const auto index_size = index_length / sizeof(RecordDataIndexSection::Entry);
    auto& data_sections   = record_->data_sections;
    data_sections.resize(index_size);
    if (in.read(data_sections.data(), index_length) != index_length) {
        throw InvalidRecord("Unexpected EOF reached");
    }
    auto index_end = atlas::io::read_struct<RecordDataIndexSection::End>(in);
    if (not index_end.valid()) {
        throw InvalidRecord("Data index section is not valid. Invalid section end marker: [" + index_end.str() + "]");
    }
    for (auto& data_section : data_sections) {
        data_section.offset += rbegin;
        if (data_section.length < sizeof(RecordDataSection::Begin) + sizeof(RecordDataSection::End)) {
            throw InvalidRecord("Unexpected data section length: [" + std::to_string(data_section.length) + "]");
        }
    }

    record_->parse();

    if (read_to_end) {
        in.seek(rbegin + r.record_length - sizeof(RecordEnd));
        RecordEnd record_end;
        read_struct(in, record_end);
        if (not record_end.valid()) {
            throw InvalidRecord("RecordEnd has unexpected format: [" + record_end.str() + "]");
        }
        if (in.position() != rbegin + r.record_length) {
            throw InvalidRecord("RecordEnd has unexpected format");
        }
    }

    return *this;
}

//---------------------------------------------------------------------------------------------------------------------

void ParsedRecord::parse() {
    for (auto& key : keys) {
        auto& item           = items.at(key);
        item.record.version_ = head.version;
        item.record.created_ = head.time;
        item.data.section(item.getInt("data.section", 0));
        item.data.endian(head.endian());
        item.data.compression(item.getString("data.compression.type", "none"));
        if (item.data.section()) {
            auto& data_section = data_sections.at(size_t(item.data.section() - 1));
            item.data.checksum(data_section.checksum);
            item.data.compressed_size(data_section.length - sizeof(RecordDataSection::Begin) -
                                      sizeof(RecordDataSection::End));
            if (item.data.compressed()) {
                item.data.size(atlas::io::uncompressed_size(item));
            }
            else {
                item.data.size(item.data.compressed_size());
            }
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
