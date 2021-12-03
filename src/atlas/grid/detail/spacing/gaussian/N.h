/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Nov 2014

#pragma once

#include <cstddef>
#include <vector>

#include "atlas/util/Factory.h"

namespace atlas {
namespace grid {
namespace spacing {
namespace gaussian {

class GaussianLatitudes {
public:
    /// @pre nlats has enough allocated memory to store the latitudes
    /// @param size of lats vector
    void assign(double lats[], const size_t size) const;

    /// @post resizes the vector to the number of latitutes
    void assign(std::vector<double>& lats) const;

    size_t N() const { return lats_.size(); }

protected:
    std::vector<double> lats_;
};

class GaussianLatitudesFactory : public util::Factory<GaussianLatitudesFactory> {
public:
    static std::string className() { return "GaussianLatitudesFactory"; }
    static const GaussianLatitudes* build(const std::string& builder) { return get(builder)->make(); }
    using Factory::Factory;

private:
    virtual const GaussianLatitudes* make() const = 0;
};

template <typename T>
class GaussianLatitudesBuilder : public GaussianLatitudesFactory {
public:
    using GaussianLatitudesFactory::GaussianLatitudesFactory;

private:
    virtual const GaussianLatitudes* make() const override { return new T(); }
};

#define DECLARE_GAUSSIAN_LATITUDES(NUMBER)       \
    class N##NUMBER : public GaussianLatitudes { \
    public:                                      \
        N##NUMBER();                             \
    };

#define LIST(...) __VA_ARGS__
#define DEFINE_GAUSSIAN_LATITUDES(NUMBER, LATS)                            \
    N##NUMBER::N##NUMBER() {                                               \
        size_t N     = NUMBER;                                             \
        double lat[] = {LATS};                                             \
        lats_.assign(lat, lat + N);                                        \
    }                                                                      \
    namespace {                                                            \
    static GaussianLatitudesBuilder<N##NUMBER> builder_N##NUMBER(#NUMBER); \
    }


DECLARE_GAUSSIAN_LATITUDES(16);
DECLARE_GAUSSIAN_LATITUDES(24);
DECLARE_GAUSSIAN_LATITUDES(32);
DECLARE_GAUSSIAN_LATITUDES(48);
DECLARE_GAUSSIAN_LATITUDES(64);
DECLARE_GAUSSIAN_LATITUDES(80);
DECLARE_GAUSSIAN_LATITUDES(96);
DECLARE_GAUSSIAN_LATITUDES(128);
DECLARE_GAUSSIAN_LATITUDES(160);
DECLARE_GAUSSIAN_LATITUDES(200);
DECLARE_GAUSSIAN_LATITUDES(256);
DECLARE_GAUSSIAN_LATITUDES(320);
DECLARE_GAUSSIAN_LATITUDES(400);
DECLARE_GAUSSIAN_LATITUDES(512);
DECLARE_GAUSSIAN_LATITUDES(576);
DECLARE_GAUSSIAN_LATITUDES(640);
DECLARE_GAUSSIAN_LATITUDES(800);
DECLARE_GAUSSIAN_LATITUDES(1024);
DECLARE_GAUSSIAN_LATITUDES(1280);
DECLARE_GAUSSIAN_LATITUDES(1600);
DECLARE_GAUSSIAN_LATITUDES(2000);
DECLARE_GAUSSIAN_LATITUDES(4000);
DECLARE_GAUSSIAN_LATITUDES(8000);

#undef DECLARE_GAUSSIAN_LATITUDES

}  // namespace gaussian
}  // namespace spacing
}  // namespace grid
}  // namespace atlas
