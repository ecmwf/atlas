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

#include <array>
#include <iosfwd>

#include "atlas/domain/detail/RectangularDomain.h"

namespace atlas {
class RectangularDomain;
class ZonalBandDomain;
}  // namespace atlas

namespace atlas {
namespace domain {

class ZonalBandDomain : public atlas::domain::RectangularLonLatDomain {
protected:
    static constexpr char units_[] = "degrees";

public:
    static bool is_global(const Interval& y);

public:
    // constructor
    ZonalBandDomain(const eckit::Parametrisation&);
    ZonalBandDomain(const Interval&);
    ZonalBandDomain(const Interval&, const double west);

    static std::string static_type() { return "zonal_band"; }
    virtual std::string type() const override { return static_type(); }

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const override;

    virtual bool global() const override { return global_; }
    virtual bool empty() const override { return (ymin() == ymax()); }

    virtual Spec spec() const override;

    virtual void print(std::ostream&) const override;

    virtual void hash(eckit::Hash&) const override;

    virtual std::string units() const override { return units_; }

    /// Check if grid includes the North pole
    virtual bool containsNorthPole() const override;

    /// Check if grid includes the South pole
    virtual bool containsSouthPole() const override;

protected:
    friend class ::atlas::RectangularDomain;

private:
    bool global_;
    double ymin_tol_;
    double ymax_tol_;
};

}  // namespace domain
}  // namespace atlas
