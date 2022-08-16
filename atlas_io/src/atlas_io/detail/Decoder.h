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

#include <memory>

#include "atlas_io/Data.h"
#include "atlas_io/Metadata.h"
#include "atlas_io/detail/TypeTraits.h"

namespace atlas {
namespace io {

class Decoder {
public:
    template <typename T, enable_if_decodable_t<T> = 0>
    explicit Decoder(T& value): self_(new DecodableItem<T>(value)) {}

    friend void decode(const atlas::io::Metadata& metadata, const atlas::io::Data& data, Decoder&);
    friend void decode(const atlas::io::Metadata& metadata, const atlas::io::Data& data, Decoder&&);

private:
    struct Decodable {
        virtual ~Decodable()                                                     = default;
        virtual void decode_(const atlas::io::Metadata&, const atlas::io::Data&) = 0;
    };

    template <typename T>
    struct DecodableItem : Decodable {
        explicit DecodableItem(T& value): data_(value) {}

        void decode_(const atlas::io::Metadata& metadata, const atlas::io::Data& encoded) override {
            decode(metadata, encoded, data_);
        }

        T& data_;
    };

    std::shared_ptr<Decodable> self_;
};

void decode(const Metadata&, const Data&, Decoder&);
void decode(const Metadata&, const Data&, Decoder&&);

}  // namespace io
}  // namespace atlas
