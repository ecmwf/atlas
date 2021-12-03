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

#include "atlas/domain/detail/ZonalBandDomain.h"

namespace atlas {
namespace domain {

class GlobalDomain : public ZonalBandDomain {
public:
    GlobalDomain();
    GlobalDomain(const double west);
    GlobalDomain(const eckit::Parametrisation& p);

    static std::string static_type() { return "global"; }
    virtual std::string type() const override { return static_type(); }

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const override { return true; }

    // Domain properties
    virtual bool global() const override { return true; }
    virtual bool empty() const override { return false; }

    virtual Spec spec() const override;

    virtual void print(std::ostream&) const override;

    virtual void hash(eckit::Hash&) const override;

    /// Check if grid includes the North pole
    virtual bool containsNorthPole() const override { return true; }

    /// Check if grid includes the South pole
    virtual bool containsSouthPole() const override { return true; }

private:
    friend class ::atlas::RectangularDomain;
};

}  // namespace domain
}  // namespace atlas
