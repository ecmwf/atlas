/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/function/SphericalHarmonic.h"

#include <atomic>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <unordered_map>
#include <utility>

#include "atlas/runtime/Exception.h"
#include "atlas/util/Constants.h"

namespace atlas::util::function {

namespace {
static double factorial(double v) {
    if (v == 0) {
        return 1;
    }
    double result = v;
    while (--v > 0) {
        result *= v;
    }
    return result;
}

static double double_factorial(double x) {
    if (x == 0 || x == -1) {
        return 1;
    }

    double result = x;
    while ((x -= 2) > 0) {
        result *= x;
    }
    return result;
}

// Associated Legendre Polynomial
static double P(const int n, const int m, const double x) {
    // No recursive calculation needed
    if (n == m) {
        return (std::pow(-1.0, m) * double_factorial(2 * m - 1) * std::pow(std::sqrt(1. - x * x), m));
    }

    if (n == m + 1) {
        return x * (2 * m + 1) * P(m, m, x);
    }

    // Formula 1
    return (x * (2 * n - 1) * P(n - 1, m, x) - (n + m - 1) * P(n - 2, m, x)) / (n - m);
}

static double K(const int n, const int m) {
    //When m is less than 0, multiply - 1 to pass in
    return std::sqrt(((2 * n + 1) * factorial(n - m)) / (4 * M_PI * factorial(n + m)));
}

template <typename Key, typename Compute>
static double cached(const Key& key, const double argument, Compute&& compute) {
    struct LastValue {
        bool initialized{};
        Key key{};
        double argument{};
        double value{};
    };
    static thread_local LastValue last;

    if (last.initialized && last.key == key && last.argument == argument) {
        return last.value;
    }

    const double value = compute();
    last = {true, key, argument, value};
    return last.value;
}

static double P_cos_COLAT(const int n, const int abs_m, const double lat) {
    return cached(std::make_pair(n, abs_m), lat, [&] {
        const double colat_rad = (90. - lat) * Constants::degreesToRadians();
        return P(n, abs_m, std::cos(colat_rad));
    });
}

}  // namespace function

struct SphericalHarmonic::Cache {
    explicit Cache(const std::uint64_t id): id(id) {}

    std::uint64_t id;
    std::unordered_map<double, double> values;
    std::mutex mutex;
};

double SphericalHarmonic::P_cos_COLAT(const double lat) const {
    struct LastValue {
        std::uint64_t cache_id{};
        double lat{};
        double value{};
    };
    static thread_local LastValue last;

    if (last.cache_id == cache_->id && last.lat == lat) {
        return last.value;
    }

    std::lock_guard<std::mutex> lock(cache_->mutex);
    const auto cached_value = cache_->values.find(lat);
    if (cached_value != cache_->values.end()) {
        last = {cache_->id, lat, cached_value->second};
        return last.value;
    }

    const double colat_rad = (90. - lat) * Constants::degreesToRadians();
    const double value     = P(n_, abs_m_, std::cos(colat_rad));
    cache_->values.emplace(lat, value);
    last = {cache_->id, lat, value};
    return last.value;
}

double spherical_harmonic(int n, int m, double lon, double lat) {
    const int abs_m = std::abs(m);

    ATLAS_ASSERT(n >= abs_m);

    if (m == 0) {
        return K(n, 0) * P_cos_COLAT(n, 0, lat);
    }
    auto K_SQRT2 = [](int n, int abs_m) {
        struct LastValue {
            int n{-1}; // Always initialize to invalid value to avoid false cache hits
            int abs_m{};
            double value{};
        };
        static thread_local LastValue last;
        if (not (last.n == n && last.abs_m == abs_m)) {
            last = {n, abs_m, K(n, abs_m) * M_SQRT2};
        }
        return last.value;
    };
    lon *= Constants::degreesToRadians();

    if (m > 0) {
        return (K_SQRT2(n, abs_m) * std::cos(abs_m * lon) * P_cos_COLAT(n, abs_m, lat));
    }
    else { // m < 0
        // When m is less than 0, multiply - 1 in advance and send it to K
        return (K_SQRT2(n, abs_m) * std::sin(abs_m * lon) * P_cos_COLAT(n, abs_m, lat));
    }
}

SphericalHarmonic::SphericalHarmonic(const int n, const int m): n_(n), m_(m), abs_m_(std::abs(m)) {
    ATLAS_ASSERT(n_ >= abs_m_);

    static std::atomic<std::uint64_t> next_cache_id{1};
    cache_ = std::make_shared<Cache>(next_cache_id.fetch_add(1, std::memory_order_relaxed));
    Knm_  = K(n_, abs_m_) * (m_ == 0 ? 1. : M_SQRT2);
}

double SphericalHarmonic::operator()(double lon, double lat) const {
    if (m_ == 0) {
        return Knm_ * P_cos_COLAT(lat);
    }

    lon *= Constants::degreesToRadians();
    if (m_ > 0) {
        return Knm_ * std::cos(m_ * lon) * P_cos_COLAT(lat);
    }
    return Knm_ * std::sin(-m_ * lon) * P_cos_COLAT(lat);
}


}  // namespace atlas::util::function
