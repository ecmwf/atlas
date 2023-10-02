/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/grid/Unstructured.h"

#include <initializer_list>
#include <iomanip>
#include <limits>
#include <memory>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/array/ArrayView.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/option.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/NormaliseLongitude.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

namespace {
static GridFactoryBuilder<Unstructured> __register_Unstructured(Unstructured::static_type());
}

namespace {
class Normalise {
public:
    Normalise(const RectangularDomain& domain):
        degrees_(domain.units() == "degrees"), normalise_(domain.xmin(), domain.xmax()) {}

    double operator()(double x) const {
        if (degrees_) {
            x = normalise_(x);
        }
        return x;
    }

private:
    const bool degrees_;
    util::NormaliseLongitude normalise_;
};
}  // namespace


Unstructured::Unstructured(const Grid& grid, Domain domain): Grid() {
    domain_ = domain;
    points_.reset(new std::vector<PointXY>);
    points_->reserve(grid.size());
    if (not domain_) {
        domain_ = GlobalDomain();
    }
    atlas::grid::IteratorXY it(grid.xy_begin());
    PointXY p;
    if (RectangularDomain(domain_)) {
        auto normalise = Normalise(RectangularDomain(domain_));
        while (it.next(p)) {
            p.x() = normalise(p.x());
            if (domain_.contains(p)) {
                points_->emplace_back(p);
            }
        }
    }
    else if (ZonalBandDomain(domain_)) {
        while (it.next(p)) {
            if (domain_.contains(p)) {
                points_->emplace_back(p);
            }
        }
    }
    else {
        while (it.next(p)) {
            points_->emplace_back(p);
        }
    }
    points_->shrink_to_fit();
}

Unstructured::Unstructured(const util::Config& config): Grid() {
    util::Config config_domain;
    if (not config.get("domain", config_domain)) {
        config_domain.set("type", "global");
    }
    domain_ = Domain(config_domain);
    std::vector<double> xy;
    if (config.get("xy", xy)) {
        const size_t N = xy.size() / 2;
        points_.reset(new std::vector<PointXY>);
        points_->reserve(N);
        for (size_t n = 0; n < N; ++n) {
            points_->emplace_back(PointXY{xy[2 * n], xy[2 * n + 1]});
        }
    }
    else {
        std::vector<double> x;
        std::vector<double> y;
        if (not config.get("x", x)) {
            throw_Exception("x missing from configuration");
        }
        if (not config.get("y", y)) {
            throw_Exception("y missing from configuration");
        }
        ATLAS_ASSERT(x.size() == y.size());
        points_.reset(new std::vector<PointXY>);
        points_->reserve(x.size());
        for (size_t n = 0; n < x.size(); ++n) {
            points_->emplace_back(PointXY{x[n], y[n]});
        }
    }
}

Unstructured::Unstructured(std::vector<PointXY>* pts): Grid(), points_(pts) {
    domain_ = GlobalDomain();
}

Unstructured::Unstructured(std::vector<PointXY>&& pts): Grid(), points_(new std::vector<PointXY>(std::move(pts))) {
    domain_ = GlobalDomain();
}

Unstructured::Unstructured(const std::vector<PointXY>& pts): Grid(), points_(new std::vector<PointXY>(pts)) {
    domain_ = GlobalDomain();
}

Unstructured::Unstructured(std::initializer_list<PointXY> initializer_list):
    Grid(), points_(new std::vector<PointXY>(initializer_list)) {
    domain_ = GlobalDomain();
}

Unstructured::Unstructured(size_t N, const double x[], const double y[], size_t xstride, size_t ystride):
 Grid(), points_(new std::vector<PointXY>(N)) {
    util::Config config_domain;
    config_domain.set("type", "global");
    domain_ = Domain(config_domain);

    std::vector<PointXY>& p = *points_;
    const idx_t npts        = static_cast<idx_t>(p.size());

    for (idx_t n = 0; n < npts; ++n) {
        p[n].assign(x[n*xstride], y[n*ystride]);
    }
}

Unstructured::Unstructured(size_t N, const double xy[]):
    Unstructured(N,xy,xy+1,2,2) {}

Unstructured::~Unstructured() = default;

Grid::uid_t Unstructured::name() const {
    if (shortName_.empty()) {
        std::ostringstream s;
        s << "unstructured." << Grid::hash().substr(0, 7);
        shortName_ = s.str();
    }
    return shortName_;
}

void Unstructured::hash(eckit::Hash& h) const {
    ATLAS_ASSERT(points_ != nullptr);

    const std::vector<PointXY>& pts = *points_;
    h.add(&pts[0], sizeof(PointXY) * pts.size());

    for (idx_t i = 0, N = static_cast<idx_t>(pts.size()); i < N; i++) {
        const PointXY& p = pts[i];
        h << p.x() << p.y();
    }

    projection().hash(h);
}

RectangularLonLatDomain Unstructured::lonlatBoundingBox() const {
    return projection_ ? projection_.lonlatBoundingBox(domain_) : domain_;
}

idx_t Unstructured::size() const {
    ATLAS_ASSERT(points_ != nullptr);
    return static_cast<idx_t>(points_->size());
}

Grid::Spec Unstructured::spec() const {
    if (cached_spec_) {
        return *cached_spec_;
    }

    cached_spec_.reset(new Grid::Spec);

    cached_spec_->set("type", static_type());

    cached_spec_->set("domain", domain().spec());
    cached_spec_->set("projection", projection().spec());

    auto it = xy_begin();
    std::vector<double> coords(2 * size());
    idx_t c(0);
    PointXY xy;
    while (it->next(xy)) {
        coords[c++] = xy.x();
        coords[c++] = xy.y();
    }

    cached_spec_->set("xy", coords);

    return *cached_spec_;
}

Grid::Config Unstructured::meshgenerator() const {
    return Config("type", "delaunay");
}
Grid::Config Unstructured::partitioner() const {
    return Config("type", "equal_regions");
}

void Unstructured::print(std::ostream& os) const {
    os << "Unstructured(Npts:" << size() << ")";
}

namespace {  // anonymous

static class unstructured : public GridBuilder {
    using Implementation = atlas::Grid::Implementation;
    using Config         = Grid::Config;

public:
    unstructured(): GridBuilder("unstructured") {}

    void print(std::ostream& os) const override {
        os << std::left << std::setw(20) << " "
           << "Unstructured grid";
    }

    const Implementation* create(const std::string& /* name */, const Config&) const override {
        throw_NotImplemented("Cannot create unstructured grid from name", Here());
    }

    const Implementation* create(const Config& config) const override { return new Unstructured(config); }

} unstructured_;

}  // anonymous namespace

extern "C" {
const Unstructured* atlas__grid__Unstructured__points(const double xy[], int shapef[], int stridesf[]) {
    size_t nb_points = shapef[1];
    ATLAS_ASSERT(shapef[0] == 2);
    size_t stride_n = stridesf[1];
    size_t stride_v = stridesf[0];
    std::vector<PointXY> points;
    points.reserve(nb_points);
    for (size_t n = 0; n < nb_points; ++n) {
        points.emplace_back(PointXY{xy[n * stride_n + 0], xy[n * stride_n + 1 * stride_v]});
    }
    return new Unstructured(std::move(points));
}

const Unstructured* atlas__grid__Unstructured__config(util::Config* conf) {
    ATLAS_ASSERT(conf != nullptr);
    const Unstructured* grid = dynamic_cast<const Unstructured*>(Grid::create(*conf));
    ATLAS_ASSERT(grid != nullptr);
    return grid;
}
}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
