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

#include <map>
#include <string>
#include <vector>


#include "atlas_io/FileStream.h"
#include "atlas_io/RecordItem.h"
#include "atlas_io/Stream.h"
#include "atlas_io/detail/Encoder.h"
#include "atlas_io/detail/NoConfig.h"
#include "atlas_io/detail/Reference.h"
#include "atlas_io/detail/TypeTraits.h"

#include "atlas_io/types/array/ArrayReference.h"

#include "atlas_io/detail/Defaults.h"


namespace atlas {
namespace io {

template <typename Interpreted, typename T>
Interpreted interprete(T& in) {
    Interpreted interpreted;
    interprete(in, interpreted);
    return interpreted;
}

//---------------------------------------------------------------------------------------------------------------------

/// @class RecordWriter
/// @brief Write record
class RecordWriter {
public:
    using Key = std::string;

public:
    /// @brief Set compression
    void compression(const std::string&);

    /// @brief Set compression off or to default
    void compression(bool);

    /// @brief Set checksum off or to default
    void checksum(bool);

    // -- set( Key, Value ) where Value can be a variety of things

    /// @brief Add link to other record item (RecordItem::URI)
    void set(const Key&, Link&&, const eckit::Configuration& = NoConfig());

    /// @brief Add item to record
    void set(const Key&, Encoder&&, const eckit::Configuration& = NoConfig());

    /// @brief Add item to record
    template <typename Value, enable_if_rvalue_t<Value&&> = 0>
    void set(const Key& key, Value&& value, const eckit::Configuration& config = NoConfig()) {
        set(key, Encoder{std::move(value)}, config);
    }

    /// @brief Add item to record
    template <typename Value>
    void set(const Key& key, const Reference<Value>& value, const eckit::Configuration& config = NoConfig()) {
        set(key, std::move(value), config);
    }

    /// @brief Add item to record
    template <typename Value, enable_if_interpretable_t<Value, ArrayReference> = 0>
    void set(const Key& key, const Value& value, const eckit::Configuration& config = NoConfig()) {
        set(key, RecordItem(interprete<ArrayReference>(value)), config);
    }

    /// @brief Add item to record
    template <typename Value, disable_if_interpretable_t<Value, ArrayReference> = 0>
    void set(const Key& key, const Value& value, const eckit::Configuration& config = NoConfig()) {
        set(key, Encoder{value}, config);
    }

    /// @brief Write new record to path
    size_t write(const eckit::PathName&, Mode = Mode::write) const;

    /// @brief Write new record to a DataHandle
    /// @pre The DataHandle must be opened for Write access.
    size_t write(eckit::DataHandle&) const;

    /// @brief Write new record to a Stream
    /// @pre The Stream must be opened for Write access.
    size_t write(Stream) const;

    /// @brief estimate maximum size of record
    ///
    /// This could be useful to write a record to a fixed size MemoryHandle
    ///
    /// @note Without compression this matches exactly the required record size.
    /// With compression active, the data sizes are assumed to be 120% of uncompressed sizes for robustness,
    /// which may seem contradictory.
    size_t estimateMaximumSize() const;

private:
    std::vector<std::string> keys_;
    std::map<std::string, Encoder> encoders_;
    std::map<std::string, DataInfo> info_;

    std::string compression_{defaults::compression_algorithm()};
    int do_checksum_{defaults::checksum_write()};
    int nb_data_sections_{0};

    std::string metadata() const;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
