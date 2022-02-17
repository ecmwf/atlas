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

#include <initializer_list>
#include <iostream>
#include <memory>
#include <string>

#include "atlas/runtime/Exception.h"
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

namespace atlas {
namespace projection {
namespace detail {

class ProjectionImpl : public util::Object {
public:
    using Spec = atlas::util::Config;
    class Jacobian : public std::array<std::array<double, 2>, 2> {
    public:
        using std::array<std::array<double, 2>, 2>::array;

        static Jacobian identity() { return Jacobian{1., 0., 0., 1.}; }

        Jacobian(double j00, double j01, double j10, double j11):
            array{std::array<double, 2>{j00, j01}, std::array<double, 2>{j10, j11}} {}

        Jacobian(std::initializer_list<std::initializer_list<double>> list):
            Jacobian{*(list.begin()->begin()), *(list.begin()->begin() + 1), *((list.begin() + 1)->begin()),
                     *((list.begin() + 1)->begin() + 1)} {}

        Jacobian operator-(const Jacobian& jac) const {
            return Jacobian{(*this)[0][0] - jac[0][0], (*this)[0][1] - jac[0][1], (*this)[1][0] - jac[1][0],
                            (*this)[1][1] - jac[1][1]};
        }

        Jacobian operator+(const Jacobian& jac) const {
            return Jacobian{(*this)[0][0] + jac[0][0], (*this)[0][1] + jac[0][1], (*this)[1][0] + jac[1][0],
                            (*this)[1][1] + jac[1][1]};
        }

        Jacobian operator*(double a) const {
            return Jacobian{(*this)[0][0] * a, (*this)[0][1] * a, (*this)[1][0] * a, (*this)[1][1] * a};
        }

        Jacobian operator*(const Jacobian& jac) const {
            return Jacobian{(*this)[0][0] * jac[0][0] + (*this)[0][1] * jac[1][0],
                            (*this)[0][0] * jac[0][1] + (*this)[0][1] * jac[1][1],
                            (*this)[1][0] * jac[0][0] + (*this)[1][1] * jac[1][0],
                            (*this)[1][0] * jac[0][1] + (*this)[1][1] * jac[1][1]};
        }

        Point2 operator*(const Point2& x) const {
            return Point2{(*this)[0][0] * x[0] + (*this)[0][1] * x[1], (*this)[1][0] * x[0] + (*this)[1][1] * x[1]};
        }

        double norm() const {
            return std::sqrt((*this)[0][0] * (*this)[0][0] + (*this)[0][1] * (*this)[0][1] +
                             (*this)[1][0] * (*this)[1][0] + (*this)[1][1] * (*this)[1][1]);
        }

        double determinant() const { return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0]; }

        Jacobian inverse() const {
            return Jacobian{(*this)[1][1], -(*this)[0][1], -(*this)[1][0], (*this)[0][0]} * (1. / determinant());
        }

        Jacobian transpose() const { return Jacobian{(*this)[1][1], (*this)[0][1], (*this)[1][0], (*this)[0][0]}; }

        double dx_dlon() const { return (*this)[JDX][JDLON]; }
        double dy_dlon() const { return (*this)[JDY][JDLON]; }
        double dx_dlat() const { return (*this)[JDX][JDLAT]; }
        double dy_dlat() const { return (*this)[JDY][JDLAT]; }

        double dlon_dx() const { return (*this)[JDLON][JDX]; }
        double dlon_dy() const { return (*this)[JDLON][JDY]; }
        double dlat_dx() const { return (*this)[JDLAT][JDX]; }
        double dlat_dy() const { return (*this)[JDLAT][JDY]; }

        friend std::ostream& operator<<(std::ostream& os, const Jacobian& jac) {
            os << jac[0][0] << " " << jac[0][1] << "\n" << jac[1][0] << " " << jac[1][1];
            return os;
        }

    private:
        enum
        {
            JDX = 0,
            JDY = 1
        };
        enum
        {
            JDLON = 0,
            JDLAT = 1
        };
    };


public:
    static const ProjectionImpl* create(const eckit::Parametrisation& p);
    static const ProjectionImpl* create(const std::string& type, const eckit::Parametrisation& p);

    ProjectionImpl()          = default;
    virtual ~ProjectionImpl() = default;  // destructor should be virtual

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
        const double normH_;
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

    void rotate(double*) const { /* do nothing */
    }
    void unrotate(double*) const { /* do nothing */
    }

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

}  // namespace detail
}  // namespace projection
}  // namespace atlas
