/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <vector>

#include "atlas/util/Checksum.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {
namespace {

template <typename Value>
util::checksum_t streamed_checksum(const std::vector<Value>& values, std::size_t split1, std::size_t split2) {
    util::checksum_t checksum{};
    util::checksum_reset(checksum);

    util::checksum_update(checksum, values.data(), split1);
    util::checksum_update(checksum, values.data() + split1, split2);
    util::checksum_update(checksum, values.data() + split1 + split2, values.size() - split1 - split2);

    return util::checksum_digest(checksum);
}

template <typename Value>
void expect_streaming_matches_oneshot(const std::vector<Value>& values, std::size_t split1, std::size_t split2) {
    EXPECT_EQ(util::checksum(values.data(), values.size()), streamed_checksum(values, split1, split2));
}

template <typename Value>
void expect_zero_length_is_supported() {
    std::vector<Value> values;

    util::checksum_t checksum{};
    util::checksum_reset(checksum);
    util::checksum_update(checksum, values.data(), values.size());

    EXPECT_EQ(util::checksum(values.data(), values.size()), 0ul);
    EXPECT_EQ(util::checksum_digest(checksum), 0ul);
}

template <typename Value>
void expect_strided_matches_selected(const std::vector<Value>& values, std::size_t count, std::size_t stride) {
    std::vector<Value> selected(count);
    for (std::size_t i = 0; i < count; ++i) {
        selected[i] = values[i * stride];
    }

    util::checksum_t checksum{};
    util::checksum_reset(checksum);
    util::checksum_update(checksum, values.data(), count, stride);

    EXPECT_EQ(util::checksum_digest(checksum), util::checksum(selected.data(), selected.size()));
}

template <typename Value>
void expect_zero_length_strided_is_supported(const std::vector<Value>& values, std::size_t stride) {
    util::checksum_t checksum{};
    util::checksum_reset(checksum);
    util::checksum_update(checksum, values.data(), 0, stride);
    EXPECT_EQ(util::checksum_digest(checksum), 0ul);
}
}  // namespace

CASE("checksum streaming matches one-shot for int") {
    std::vector<int> values(257);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<int>(i * 17 + 3);
    }

    expect_streaming_matches_oneshot(values, 13, 101);
}

CASE("checksum streaming matches one-shot for long") {
    std::vector<long> values(129);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<long>(i * 101 + 7);
    }

    expect_streaming_matches_oneshot(values, 17, 53);
}

CASE("checksum streaming matches one-shot for float") {
    std::vector<float> values(73);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = 0.25f * static_cast<float>(i + 1) * static_cast<float>(i + 5);
    }

    expect_streaming_matches_oneshot(values, 11, 29);
}

CASE("checksum streaming matches one-shot for double") {
    std::vector<double> values(37);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = 0.5 * static_cast<double>(i + 1) * static_cast<double>(i + 3);
    }

    expect_streaming_matches_oneshot(values, 5, 11);
}

CASE("checksum streaming matches one-shot for checksum_t") {
    std::vector<util::checksum_t> values(19);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<util::checksum_t>(i * 1009u + 11u);
    }

    expect_streaming_matches_oneshot(values, 4, 7);
}

CASE("checksum supports zero-length int input") {
    expect_zero_length_is_supported<int>();
}

CASE("checksum supports zero-length long input") {
    expect_zero_length_is_supported<long>();
}

CASE("checksum supports zero-length float input") {
    expect_zero_length_is_supported<float>();
}

CASE("checksum supports zero-length double input") {
    expect_zero_length_is_supported<double>();
}

CASE("checksum supports zero-length checksum_t input") {
    expect_zero_length_is_supported<util::checksum_t>();
}

CASE("checksum strided update matches selected values for int") {
    std::vector<int> values(300);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<int>(i * 13 + 5);
    }

    expect_strided_matches_selected(values, 73, 3);
    expect_strided_matches_selected(values, 100, 1);
    expect_zero_length_strided_is_supported(values, 7);
}

CASE("checksum strided update matches selected values for long") {
    std::vector<long> values(260);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<long>(i * 71 + 17);
    }

    expect_strided_matches_selected(values, 52, 5);
    expect_strided_matches_selected(values, 87, 1);
    expect_zero_length_strided_is_supported(values, 9);
}

CASE("checksum strided update matches selected values for float") {
    std::vector<float> values(210);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = 0.125f * static_cast<float>((i + 2) * (i + 7));
    }

    expect_strided_matches_selected(values, 42, 4);
    expect_strided_matches_selected(values, 64, 1);
    expect_zero_length_strided_is_supported(values, 11);
}

CASE("checksum strided update matches selected values for double") {
    std::vector<double> values(205);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = 0.0625 * static_cast<double>((i + 3) * (i + 9));
    }

    expect_strided_matches_selected(values, 41, 5);
    expect_strided_matches_selected(values, 70, 1);
    expect_zero_length_strided_is_supported(values, 13);
}

CASE("checksum strided update matches selected values for checksum_t") {
    std::vector<util::checksum_t> values(190);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<util::checksum_t>(i * 4099u + 29u);
    }

    expect_strided_matches_selected(values, 38, 5);
    expect_strided_matches_selected(values, 55, 1);
    expect_zero_length_strided_is_supported(values, 15);
}

CASE("checksum reset clears state") {
    std::vector<int> values(32);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<int>(i * 9 + 1);
    }

    util::checksum_t checksum{};
    util::checksum_reset(checksum);
    util::checksum_update(checksum, values.data(), 10);
    util::checksum_reset(checksum);
    util::checksum_update(checksum, values.data(), values.size());

    EXPECT_EQ(util::checksum_digest(checksum), util::checksum(values.data(), values.size()));
}

CASE ("test_vector") {
    std::vector<int> int_vector = {1, 2};
    EXPECT_EQ( util::checksum(int_vector.data(), int_vector.size()), 0x1003);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}