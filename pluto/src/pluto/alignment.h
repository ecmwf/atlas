/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstddef>
#include <cstdint>

namespace pluto {

// --------------------------------------------------------------------------------------------------------

constexpr std::size_t default_alignment() {
    return 256;  // This is device-arch dependent. hicMalloc & co don't have alignment argument
}

// --------------------------------------------------------------------------------------------------------

/**
 * @brief Returns whether or not `alignment` is a valid memory alignment.
 *
 * @param[in] alignment to check
 *
 * @return True if the alignment is valid, false otherwise.
 */
[[nodiscard]] constexpr bool is_supported_alignment(std::size_t alignment) noexcept {
    constexpr auto is_pow2 = [](std::size_t value) { return (value != 0U) && ((value & (value - 1)) == 0U); };
    return is_pow2(alignment);
}

/**
 * @brief Align up to nearest multiple of specified power of 2
 *
 * @param[in] value value to align
 * @param[in] alignment amount, in bytes, must be a power of 2
 *
 * @return the aligned value
 */
[[nodiscard]] constexpr std::size_t align_up(std::size_t value, std::size_t alignment) noexcept {
    return (value + (alignment - 1)) & ~(alignment - 1);
}

/**
 * @brief Align down to the nearest multiple of specified power of 2
 *
 * @param[in] value value to align
 * @param[in] alignment amount, in bytes, must be a power of 2
 *
 * @return the aligned value
 */
[[nodiscard]] constexpr std::size_t align_down(std::size_t value, std::size_t alignment) noexcept {
    return value & ~(alignment - 1);
}

/**
 * @brief Checks whether a value is aligned to a multiple of a specified power of 2
 *
 * @param[in] value value to check for alignment
 * @param[in] alignment amount, in bytes, must be a power of 2
 *
 * @return true if aligned
 */
[[nodiscard]] constexpr bool is_aligned(std::size_t value, std::size_t alignment) noexcept {
    return value == align_down(value, alignment);
}

/**
 * @brief Checks whether the provided pointer is aligned to a specified @p alignment
 *
 * @param[in] ptr pointer to check for alignment
 * @param[in] alignment required alignment in bytes, must be a power of 2
 *
 * @return true if the pointer is aligned
 */
[[nodiscard]] inline bool is_aligned(void* ptr, std::size_t alignment) noexcept {
    // return is_aligned(reinterpret_cast<std::uintptr_t>(ptr), alignment);
    return reinterpret_cast<std::uintptr_t>(ptr) % alignment == 0;
}


// --------------------------------------------------------------------------------------------------------

}  // namespace pluto
