/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "eckit/geometry/Point2.h"
#include "eckit/log/Plural.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/domain.h"
#include "atlas/domain/detail/GlobalDomain.h"  // FIXME not included by atlas/domain.h
#include "atlas/grid.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Point.h"
#include "atlas/util/Rotation.h"
#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

double normalise(double weird, double minimum) {
    double lon = weird;

    double GLOBE = 360.;
    while (lon < minimum) {
        lon += GLOBE;
    }
    while (lon >= minimum + GLOBE) {
        lon -= GLOBE;
    }

    return lon;
}


struct BoundingBox : public std::array<double, 4> {
    BoundingBox(double n, double w, double s, double e): std::array<double, 4>({n, w, s, e}) {}
    BoundingBox(): BoundingBox(90, 0, -90, 360) {}

    double north() const { return operator[](0); }
    double south() const { return operator[](2); }
    double west() const { return operator[](1); }
    double east() const { return operator[](3); }

    operator RectangularDomain() const { return {{{west(), east()}}, {{south(), north()}}}; }

    bool contains(const double& lat, const double& lon) const {
        return (lat <= north()) && (lat >= south()) && (normalise(lon, west()) <= east());
    }

    bool contains(const BoundingBox& other) const {
        bool otherEmpty = !eckit::types::is_strictly_greater(other.north(), other.south()) ||
                          !eckit::types::is_strictly_greater(other.east(), other.west());

        if (otherEmpty) {
            return contains(other.south(), other.west());
        }

        // check for West/East range (if non-periodic), then other's corners
        if (east() - west() < other.east() - other.west() || east() < normalise(other.east(), west())) {
            return false;
        }

        return contains(other.north(), other.west()) && contains(other.north(), other.east()) &&
               contains(other.south(), other.west()) && contains(other.south(), other.east());
    }
};


struct Rotation : std::array<double, 2> {
    Rotation(double lat, double lon): std::array<double, 2>({lon, lat}) {}
    double south_pole_latitude() const { return operator[](1); }
    double south_pole_longitude() const { return normalise(operator[](0), 0.); }
    BoundingBox rotate(const BoundingBox& box) const {
        util::Config config;
        config.set("type", "rotated_lonlat");
        config.set("south_pole", std::vector<double>({south_pole_longitude(), south_pole_latitude()}));
        atlas::Projection p(config);

        RectangularDomain before({box.west(), box.east()}, {box.south(), box.north()});
        Domain after = p.lonlatBoundingBox(before);
        ATLAS_ASSERT(after);

        RectangularLonLatDomain r(after);
        ATLAS_ASSERT(r);

        return {r.north(), r.west(), r.south(), r.east()};
    }
};


using NiNj = std::array<long, 2>;


struct test_poles_t {
public:
    test_poles_t(const NiNj& ninj, const Rotation& rotation, const BoundingBox& bbox, bool includesNorthPole = false,
                 bool includesSouthPole = false):
        ninj_(ninj),
        rotation_(rotation),
        bbox_(bbox),
        includesNorthPole_(includesNorthPole),
        includesSouthPole_(includesSouthPole) {}

    const NiNj ninj_;
    const Rotation rotation_;
    const BoundingBox bbox_;
    const bool includesNorthPole_;
    const bool includesSouthPole_;

private:
    friend std::ostream& operator<<(std::ostream& out, const test_poles_t& test) {
        return (out << "test:"
                    << "\n\t"
                    << "NiNj        = " << test.ninj_ << "\n\t"
                    << "Rotation    = " << test.rotation_ << "\n\t"
                    << "BoundingBox = " << test.bbox_ << "\n\t"
                    << "includesNorthPole? " << test.includesNorthPole_ << "\n\t"
                    << "includesSouthPole? " << test.includesSouthPole_);
    }
};


//-----------------------------------------------------------------------------

CASE("MIR-282") {
    auto old = Log::info().precision(16);


    std::vector<test_poles_t> test_poles{
        {NiNj{124, 118}, Rotation(-35., 0.), BoundingBox{12, -14.5, -17.25, 16.25}},
        {NiNj{360, 181}, Rotation(-90., 0.), BoundingBox(), true, true},
        {NiNj{81, 56}, Rotation(-75., 15.), BoundingBox{75.1, -35., 20., 45.}, true},
        {NiNj{111, 86}, Rotation(-35., 15.), BoundingBox{40., -55., -45., 55.}, true},
        {NiNj{91, 76}, Rotation(-30., -15.), BoundingBox{35., -40., -40., 50.}, true},
        {NiNj{101, 81}, Rotation(-25., 0.), BoundingBox{40., -50., -40., 50.}, true},
        {NiNj{56, 61}, Rotation(-15., 45.), BoundingBox{30., -50., -30., 5.}, true},
        {NiNj{96, 91}, Rotation(0., 80.), BoundingBox{50., -65., -40., 30.}, true},

        {NiNj{178, 143}, Rotation(-40., 10.), BoundingBox{22.7, -13.6, -5.9, 21.8}},
        {NiNj{117, 79}, Rotation(-43., 10.), BoundingBox{3.4, -6.8, -4.4, 4.8}},
        {NiNj{776, 492}, Rotation(-30., 0.), BoundingBox{18.1, -37.6, -31., 39.9}},
        {NiNj{149, 105}, Rotation(-76., 14.), BoundingBox{72., -32., 20., 42.}},
        {NiNj{240, 240}, Rotation(-30., -5.), BoundingBox{9.875, -15., -20., 14.875}},

        {NiNj{13, 12}, Rotation(-15., 45.), BoundingBox{27.5, -46.5, -28, 1.5}, true},
        {NiNj{321, 370}, Rotation(-15., 45.), BoundingBox{27.5, -46.5, -28, 1.5}, true},

        {NiNj{202, 235}, Rotation(0., 130.), BoundingBox{32.75, -86.75, -37.75, -26.15}, false},
        {NiNj{409, 309}, Rotation(-35., 15.), BoundingBox{36., -51., -41., 51.}, true},

        {NiNj{171, 6}, Rotation(-35., 160.), BoundingBox{80., 30., 75., 200.}},
        {NiNj{81, 11}, Rotation(30., -30.), BoundingBox{70., 120., 60., 200.}},
        {NiNj{261, 66}, Rotation(45., -120.), BoundingBox{55., -120., -10., 140.}},

        {NiNj{4, 4}, Rotation(50., 100.), BoundingBox{10., 70., -20., 100.}},
    };


    SECTION("MIR-282: rotated_ll covering North/South poles") {
        for (auto& test : test_poles) {
            Log::info() << '\n' << test << std::endl;

            const PointLonLat southPole(test.rotation_.south_pole_longitude(), test.rotation_.south_pole_latitude());

            const util::Rotation r(southPole);

            // check bbox including poles (in the unrotated frame)
            PointLonLat NP{r.unrotate({0., 90.})};
            PointLonLat SP{r.unrotate({0., -90.})};

            bool includesNorthPole = test.bbox_.contains(NP.lat(), NP.lon());
            bool includesSouthPole = test.bbox_.contains(SP.lat(), SP.lon());

            Log::info() << "check:"
                        << "\n\t"
                        << "NP = " << NP << "\n\t"
                        << "SP = " << SP << "\n\t"
                        << "includesNorthPole? " << includesNorthPole << "\n\t"
                        << "includesSouthPole? " << includesSouthPole << std::endl;

            EXPECT(includesNorthPole == test.includesNorthPole_);
            EXPECT(includesSouthPole == test.includesSouthPole_);
        }
    }


    SECTION("MIR-282: rotated_ll contained by cropping") {
        for (auto& test : test_poles) {
            Log::info() << '\n' << test << std::endl;

            double n = test.bbox_.north();
            double s = test.bbox_.south();
            double w = test.bbox_.west();
            double e = test.bbox_.east();

            RectangularDomain dom(test.bbox_);

            using grid::LinearSpacing;
            StructuredGrid::XSpace xspace(LinearSpacing(w, e, test.ninj_[0], !ZonalBandDomain(dom)));
            StructuredGrid::YSpace yspace(LinearSpacing(n, s, test.ninj_[1]));

            util::Config config;
            config.set("type", "rotated_lonlat");
            config.set("south_pole", std::vector<double>({test.rotation_.south_pole_longitude(),
                                                          test.rotation_.south_pole_latitude()}));
            Projection rotation(config);

            Grid g = StructuredGrid(xspace, yspace, rotation, dom);
            ATLAS_ASSERT(g);

            BoundingBox crop(test.rotation_.rotate(test.bbox_));

            Log::info() << "contained by cropping"
                        << "\n\t   " << test.bbox_ << "\n\t + " << test.rotation_ << "\n\t = " << crop << std::endl;

            for (PointLonLat p : g.lonlat()) {
                EXPECT(crop.contains(p.lat(), p.lon()));
            }
            Log::info() << "checked " << eckit::Plural(g.size(), "point") << std::endl;
        }
    }


    Log::info().precision(old);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
