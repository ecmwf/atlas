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

#include <iosfwd>

#include "atlas/util/Config.h"
#include "atlas/util/Object.h"
#include "atlas/util/Point.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
class Projection;

namespace domain {

class Domain : public util::Object {
public:
    using Spec = util::Config;

public:
    static const Domain* create();  // Create a global domain

    static const Domain* create(const eckit::Parametrisation&);

    virtual std::string type() const = 0;

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const = 0;

    bool contains(const PointXY& p) const { return contains(p.x(), p.y()); }

    // Specification of grid
    virtual Spec spec() const = 0;

    /// Check if domain represents the complete globe surface
    virtual bool global() const = 0;

    /// Check if domain does not represent any area on the globe surface
    virtual bool empty() const = 0;

    // Unless the domain is global, we can never be sure about these functions
    // without knowing also the projection

    /// Check if grid includes the North pole
    virtual bool containsNorthPole() const = 0;

    /// Check if grid includes the South pole
    virtual bool containsSouthPole() const = 0;

    /// Output to stream
    virtual void print(std::ostream&) const = 0;

    friend std::ostream& operator<<(std::ostream& s, const Domain& d) {
        d.print(s);
        return s;
    }

    virtual std::string units() const = 0;

    virtual void hash(eckit::Hash&) const = 0;
};

class RectangularDomain;

extern "C" {
const Domain* atlas__Domain__ctor_config(const eckit::Parametrisation* config);
void atlas__Domain__type(const Domain* This, char*& type, int& size);
void atlas__Domain__hash(const Domain* This, char*& hash, int& size);
Domain::Spec* atlas__Domain__spec(const Domain* This);
double atlas__LonLatRectangularDomain__north(const RectangularDomain* This);
double atlas__LonLatRectangularDomain__west(const RectangularDomain* This);
double atlas__LonLatRectangularDomain__south(const RectangularDomain* This);
double atlas__LonLatRectangularDomain__east(const RectangularDomain* This);
}

}  // namespace domain
}  // namespace atlas
