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
#include <thread>
#include <vector>

#include "tests/AtlasTestEnvironment.h"

namespace atlas::test {

namespace {

double expected(int n, int m, double lon, double lat) {
    constexpr double pi      = M_PI;
    lon *= pi / 180.;
    lat *= pi / 180.;

    if (n == 0 && m == 0) {
        return 1. / std::sqrt(4. * pi);
    }
    if (n == 1 && m == 0) {
        return std::sqrt(3. / (4. * pi)) * std::sin(lat);
    }
    if (n == 1 && m == 1) {
        return -std::sqrt(3. / (4. * pi)) * std::cos(lat) * std::cos(lon);
    }
    if (n == 1 && m == -1) {
        return -std::sqrt(3. / (4. * pi)) * std::cos(lat) * std::sin(lon);
    }
    if (n == 2 && m == 1) {
        return -3. * std::sqrt(5. / (12. * pi)) * std::sin(lat) * std::cos(lat) *
               std::cos(lon);
    }
    if (n == 2 && m == -1) {
        return -3. * std::sqrt(5. / (12. * pi)) * std::sin(lat) * std::cos(lat) *
               std::sin(lon);
    }
    ATLAS_NOTIMPLEMENTED;
}

void expect_harmonic(int n, int m, double lon, double lat) {
    const double value = expected(n, m, lon, lat);
    EXPECT_APPROX_EQ(util::function::spherical_harmonic(n, m, lon, lat), value, 1.e-14);
    EXPECT_APPROX_EQ(util::function::SphericalHarmonic(n, m)(lon, lat), value, 1.e-14);
}

}  // namespace

CASE("spherical harmonic function and class match analytic values") {
    expect_harmonic(0, 0, 0., 0.);
    expect_harmonic(1, 0, 45., 30.);
    expect_harmonic(1, 1, 35., -20.);
    expect_harmonic(1, -1, 35., -20.);
    expect_harmonic(2, 1, -70., 25.);
    expect_harmonic(2, -1, -70., 25.);
}

CASE("spherical harmonic caches handle latitude reuse and pair changes") {
    const util::function::SphericalHarmonic first(1, 1);
    const util::function::SphericalHarmonic second(2, 1);

    for (const double lat : {20., -35., 20.}) {
        EXPECT_APPROX_EQ(first(40., lat), expected(1, 1, 40., lat), 1.e-14);
        EXPECT_APPROX_EQ(util::function::spherical_harmonic(1, 1, 40., lat), expected(1, 1, 40., lat), 1.e-14);
    }

    for (int iteration = 0; iteration < 3; ++iteration) {
        EXPECT_APPROX_EQ(first(40., 20.), expected(1, 1, 40., 20.), 1.e-14);
        EXPECT_APPROX_EQ(second(40., 20.), expected(2, 1, 40., 20.), 1.e-14);
        EXPECT_APPROX_EQ(first(40., 20.), expected(1, 1, 40., 20.), 1.e-14);
    }
}

CASE("spherical harmonic rejects invalid wave numbers") {
    EXPECT_THROWS(util::function::spherical_harmonic(2, 3, 0., 0.));
    EXPECT_THROWS(util::function::SphericalHarmonic(2, 3));
}

CASE("spherical harmonic caches are thread safe") {
    constexpr int thread_count = 8;
    const util::function::SphericalHarmonic first(1, 1);
    const util::function::SphericalHarmonic second(2, -1);
    std::atomic<bool> start{false};
    std::vector<int> failures(thread_count, 0);
    std::vector<std::thread> threads;
    threads.reserve(thread_count);

    for (int thread = 0; thread < thread_count; ++thread) {
        threads.emplace_back([thread, &first, &second, &start, &failures]() {
            while (!start.load(std::memory_order_acquire)) {
                std::this_thread::yield();
            }

            for (int iteration = 0; iteration < 100; ++iteration) {
                const double lon = -150. + 7. * thread + iteration;
                const double lat = -60. + 5. * thread + (iteration % 5) * 11.;
                failures[thread] += std::abs(first(lon, lat) - expected(1, 1, lon, lat)) > 1.e-14;
                failures[thread] += std::abs(second(lon, lat) - expected(2, -1, lon, lat)) > 1.e-14;
                failures[thread] +=
                    std::abs(util::function::spherical_harmonic(1, -1, lon, lat) - expected(1, -1, lon, lat)) >
                    1.e-14;
            }
        });
    }

    start.store(true, std::memory_order_release);
    for (auto& thread : threads) {
        thread.join();
    }

    for (const int failure_count : failures) {
        EXPECT_EQ(failure_count, 0);
    }
}

}  // namespace atlas::test

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}