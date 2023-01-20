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
#include <iosfwd>
#include <memory>
#include <string>

#include "atlas_io/Data.h"
#include "atlas_io/RecordItem.h"
#include "atlas_io/detail/DataInfo.h"
#include "atlas_io/detail/Link.h"
#include "atlas_io/detail/Reference.h"
#include "atlas_io/detail/TypeTraits.h"

namespace atlas {
namespace io {

class Encoder {
public:
    Encoder() = default;

    operator bool() const { return bool(self_); }

    template <typename T, enable_if_move_constructible_encodable_rvalue_t<T> = 0>
    explicit Encoder(T&& x): self_(new EncodableValue<T>(std::move(x))) {}

    Encoder(const Link& link): self_(new EncodableLink(link)) {}

    Encoder(Encoder&& other): self_(std::move(other.self_)) {}

    template <typename T, enable_if_scalar_t<T> = 0>
    explicit Encoder(const T& x): self_(new EncodableValue<T>(x)) {}


    Encoder& operator=(Encoder&& rhs) {
        self_ = std::move(rhs.self_);
        return *this;
    }

    friend size_t encode_metadata(const Encoder&, atlas::io::Metadata&);

    friend void encode_data(const Encoder&, atlas::io::Data&);

    bool encodes_data() const { return self_->encodes_data_(); }

private:
    struct Encodable {
        virtual ~Encodable()                                        = default;
        virtual size_t encode_metadata_(atlas::io::Metadata&) const = 0;
        virtual void encode_data_(atlas::io::Data&) const           = 0;
        virtual bool encodes_data_() const                          = 0;
    };

    template <typename Value>
    struct EncodableValue : Encodable {
        EncodableValue(Value&& v): value_{std::move(v)} { sfinae::encode_metadata(value_, metadata_, data_size_); }

        template <bool EnableBool = true, enable_if_scalar_t<Value, EnableBool> = 0>
        EncodableValue(const Value& v): value_{v} {
            sfinae::encode_metadata(value_, metadata_, data_size_);
        }

        size_t encode_metadata_(atlas::io::Metadata& metadata) const override {
            metadata.set(metadata_);
            return data_size_;
        }

        void encode_data_(atlas::io::Data& out) const override { sfinae::encode_data(value_, out); }
        bool encodes_data_() const override { return data_size_ > 0; }

        const Value value_;
        Metadata metadata_;
        size_t data_size_{0};
    };

    struct EncodableLink : Encodable {
        EncodableLink(const Link& link): link_(link) {}

        size_t encode_metadata_(atlas::io::Metadata& metadata) const override {
            metadata.set(atlas::io::Metadata("link", link_.uri));
            return 0;
        }

        void encode_data_(atlas::io::Data& /*out*/) const override {}
        bool encodes_data_() const override { return false; }

        Link link_;
    };


    std::unique_ptr<Encodable> self_;
};

size_t encode_metadata(const Encoder& encoder, atlas::io::Metadata& metadata);
void encode_data(const Encoder& encoder, atlas::io::Data& out);


}  // namespace io
}  // namespace atlas
