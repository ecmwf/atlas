/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "atlas/projection/Jacobian.h"
#include "atlas/util/Factory.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Object.h"
#include "atlas/util/Point.h"
#include "atlas/util/Rotation.h"

namespace eckit {
class Parametrisation;
class Hash;
}  // namespace eckit

namespace atlas {
class Domain;
class RectangularLonLatDomain;
namespace util {
class Config;
}
}  // namespace atlas

namespace atlas::projection::detail {

class ProjectionImpl : public util::Object {
public:
    using Spec     = atlas::util::Config;
    using Jacobian = projection::Jacobian;

public:
    static const ProjectionImpl* create(const eckit::Parametrisation& p);
    static const ProjectionImpl* create(const std::string& type, const eckit::Parametrisation& p);

    ProjectionImpl() = default;

    virtual std::string type() const = 0;

    virtual void xy2lonlat(double crd[]) const = 0;
    virtual void lonlat2xy(double crd[]) const = 0;

    virtual Jacobian jacobian(const PointLonLat&) const = 0;

    void xy2lonlat(Point2&) const;
    void lonlat2xy(Point2&) const;

    PointLonLat lonlat(const PointXY&) const;
    PointXY xy(const PointLonLat&) const;
    virtual PointXYZ xyz(const PointLonLat&) const;

    virtual bool strictlyRegional() const                                  = 0;
    virtual RectangularLonLatDomain lonlatBoundingBox(const Domain&) const = 0;

    virtual Spec spec() const = 0;

    virtual std::string units() const = 0;

    virtual operator bool() const { return true; }

    virtual void hash(eckit::Hash&) const = 0;

    struct BoundLonLat {
        operator RectangularLonLatDomain() const;
        void extend(PointLonLat p, PointLonLat eps);

        bool crossesDateLine(bool);
        bool includesNorthPole(bool);
        bool includesSouthPole(bool);

        bool crossesDateLine() const { return crossesDateLine_; }
        bool includesNorthPole() const { return includesNorthPole_; }
        bool includesSouthPole() const { return includesSouthPole_; }

    private:
        PointLonLat min_;
        PointLonLat max_;
        bool crossesDateLine_   = false;
        bool includesNorthPole_ = false;
        bool includesSouthPole_ = false;
        bool first_             = true;
    };

    struct Normalise {
        Normalise(const eckit::Parametrisation&);
        Normalise(double west);
        void hash(eckit::Hash&) const;
        void spec(Spec&) const;
        void operator()(double crd[]) const {
            if (normalise_) {
                crd[0] = (*normalise_)(crd[0]);
            }
        }
        operator bool() const { return bool(normalise_); }

    private:
        std::unique_ptr<util::NormaliseLongitude> normalise_;
        std::vector<double> values_;
    };

    struct Derivate {
        Derivate(const ProjectionImpl& p, PointXY A, PointXY B, double h, double refLongitude = 0.);
        virtual ~Derivate();
        virtual PointLonLat d(PointXY) const = 0;

    protected:
        const ProjectionImpl& projection_;
        const PointXY H_;
        const double invnH_;
        const double refLongitude_;
        PointLonLat xy2lonlat(const PointXY& p) const;
    };

    struct DerivateFactory : public util::Factory<DerivateFactory> {
        static std::string className() { return "DerivateFactory"; }
        static ProjectionImpl::Derivate* build(const std::string& type, const ProjectionImpl& p, PointXY A, PointXY B,
                                               double h, double refLongitude = 0.);

    protected:
        using Factory::Factory;
        virtual ~DerivateFactory();
        virtual ProjectionImpl::Derivate* make(const ProjectionImpl& p, PointXY A, PointXY B, double h,
                                               double refLongitude = 0.) = 0;
    };
};

inline void ProjectionImpl::xy2lonlat(Point2& point) const {
    xy2lonlat(point.data());
}

inline void ProjectionImpl::lonlat2xy(Point2& point) const {
    lonlat2xy(point.data());
}

inline PointLonLat ProjectionImpl::lonlat(const PointXY& xy) const {
    PointLonLat lonlat(xy);
    xy2lonlat(lonlat.data());
    return lonlat;
}

inline PointXY ProjectionImpl::xy(const PointLonLat& lonlat) const {
    PointXY xy(lonlat);
    lonlat2xy(xy.data());
    return xy;
}

//---------------------------------------------------------------------------------------------------------------------

class Rotated : public util::Rotation {
public:
    using Spec = ProjectionImpl::Spec;

    Rotated(const PointLonLat& south_pole, double rotation_angle = 0.);
    Rotated(const eckit::Parametrisation&);
    virtual ~Rotated() = default;

    static std::string classNamePrefix() { return "Rotated"; }
    static std::string typePrefix() { return "rotated_"; }

    void spec(Spec&) const;

    void hash(eckit::Hash&) const;
};

class NotRotated {
public:
    using Spec = ProjectionImpl::Spec;

    NotRotated() = default;
    NotRotated(const eckit::Parametrisation&) {}
    virtual ~NotRotated() = default;

    static std::string classNamePrefix() { return ""; }  // deliberately empty
    static std::string typePrefix() { return ""; }       // deliberately empty

    void rotate(double*) const { /* do nothing */ }
    void unrotate(double*) const { /* do nothing */ }

    bool rotated() const { return false; }

    void spec(Spec&) const {}

    void hash(eckit::Hash&) const {}
};

extern "C" {
const ProjectionImpl* atlas__Projection__ctor_config(const eckit::Parametrisation* config);
void atlas__Projection__type(const ProjectionImpl* This, char*& type, int& size);
void atlas__Projection__hash(const ProjectionImpl* This, char*& hash, int& size);
ProjectionImpl::Spec* atlas__Projection__spec(const ProjectionImpl* This);
void atlas__Projection__xy2lonlat(const ProjectionImpl* This, const double x, const double y, double& lon, double& lat);
void atlas__Projection__lonlat2xy(const ProjectionImpl* This, const double lon, const double lat, double& x, double& y);
}

}  // namespace atlas::projection::detail
