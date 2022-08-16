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
#include <memory>
#include <string>

#include "atlas_io/Data.h"
#include "atlas_io/Metadata.h"
#include "atlas_io/detail/tag.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

struct RecordItem {
public:
    struct URI {
        URI(const std::string& uri);
        URI(const std::string& _path, std::uint64_t _offset, const std::string& _key);
        std::string str() const;
        operator std::string() { return str(); }
        std::string path;
        std::uint64_t offset;
        std::string key;
    };

public:
    RecordItem() = default;

    template <typename T>
    explicit RecordItem(T&& x, tag::enable_static_assert = tag::enable_static_assert());

    template <typename T>
    explicit RecordItem(T&& x, tag::disable_static_assert);

    RecordItem(Metadata&&, Data&&);

    RecordItem(RecordItem&& other);

    const Data& data() const;

    const Metadata& metadata() const;

    void metadata(const Metadata& m);

    void data(atlas::io::Data&& d);

    bool empty() const;

    void clear();

    void decompress();

    void compress();


private:
    std::unique_ptr<Metadata> metadata_{new Metadata()};
    Data data_;
};

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
RecordItem::RecordItem(T&& x, tag::enable_static_assert) {
    encode(x, *metadata_, data_);
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
RecordItem::RecordItem(T&& x, tag::disable_static_assert) {
    encode(x, *metadata_, data_, tag::disable_static_assert());
}

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata(const RecordItem& in, atlas::io::Metadata& metadata);

//---------------------------------------------------------------------------------------------------------------------

void encode_data(const RecordItem& in, atlas::io::Data& out);

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
