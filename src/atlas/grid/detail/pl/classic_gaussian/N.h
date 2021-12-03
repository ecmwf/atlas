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
/// @date Mar 2016

#pragma once

#include <cstddef>
#include <vector>

#include "atlas/util/Factory.h"

namespace atlas {
namespace grid {
namespace detail {
namespace pl {
namespace classic_gaussian {

class PointsPerLatitude {
public:
    /// @pre nlats has enough allocated memory to store the latitudes
    /// @param size of lats vector
    void assign(long nlon[], const size_t size) const;

    /// @pre nlats has enough allocated memory to store the latitudes
    /// @param size of lats vector
    void assign(int nlon[], const size_t size) const;

    /// @post resizes the vector to the number of latitutes
    void assign(std::vector<long>& nlon) const;

    /// @post resizes the vector to the number of latitutes
    void assign(std::vector<int>& nlon) const;

    size_t N() const { return nlon_.size(); }

protected:
    std::vector<long> nlon_;
};


class PointsPerLatitudeFactory : public util::Factory<PointsPerLatitudeFactory> {
public:
    static std::string className() { return "PointsPerLatitudeFactory"; }
    static const PointsPerLatitude* build(const std::string& builder) { return get(builder)->make(); }
    using Factory::Factory;

private:
    virtual const PointsPerLatitude* make() const = 0;
};

template <typename T>
class PointsPerLatitudeBuilder : public PointsPerLatitudeFactory {
public:
    using PointsPerLatitudeFactory::PointsPerLatitudeFactory;

private:
    virtual const PointsPerLatitude* make() const override { return new T(); }
};


#define DECLARE_POINTS_PER_LATITUDE(NUMBER)      \
    class N##NUMBER : public PointsPerLatitude { \
    public:                                      \
        N##NUMBER();                             \
    };

#define LIST(...) __VA_ARGS__
#define DEFINE_POINTS_PER_LATITUDE(NUMBER, NLON)                           \
    N##NUMBER::N##NUMBER() {                                               \
        size_t N    = NUMBER;                                              \
        long nlon[] = {NLON};                                              \
        nlon_.assign(nlon, nlon + N);                                      \
    }                                                                      \
    namespace {                                                            \
    static PointsPerLatitudeBuilder<N##NUMBER> builder_N##NUMBER(#NUMBER); \
    }


DECLARE_POINTS_PER_LATITUDE(16);
DECLARE_POINTS_PER_LATITUDE(24);
DECLARE_POINTS_PER_LATITUDE(32);
DECLARE_POINTS_PER_LATITUDE(48);
DECLARE_POINTS_PER_LATITUDE(64);
DECLARE_POINTS_PER_LATITUDE(80);
DECLARE_POINTS_PER_LATITUDE(96);
DECLARE_POINTS_PER_LATITUDE(128);
DECLARE_POINTS_PER_LATITUDE(160);
DECLARE_POINTS_PER_LATITUDE(200);
DECLARE_POINTS_PER_LATITUDE(256);
DECLARE_POINTS_PER_LATITUDE(320);
DECLARE_POINTS_PER_LATITUDE(400);
DECLARE_POINTS_PER_LATITUDE(512);
DECLARE_POINTS_PER_LATITUDE(576);
DECLARE_POINTS_PER_LATITUDE(640);
DECLARE_POINTS_PER_LATITUDE(800);
DECLARE_POINTS_PER_LATITUDE(1024);
DECLARE_POINTS_PER_LATITUDE(1280);
DECLARE_POINTS_PER_LATITUDE(1600);
DECLARE_POINTS_PER_LATITUDE(2000);
DECLARE_POINTS_PER_LATITUDE(4000);
DECLARE_POINTS_PER_LATITUDE(8000);

#undef DECLARE_POINTS_PER_LATITUDE

}  // namespace classic_gaussian
}  // namespace pl
}  // namespace detail
}  // namespace grid
}  // namespace atlas
