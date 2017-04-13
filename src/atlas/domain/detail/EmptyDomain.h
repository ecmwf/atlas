#pragma once

#include "atlas/domain/detail/Domain.h"

namespace atlas {
namespace domain {

class EmptyDomain: public Domain {

public:

    EmptyDomain();
    EmptyDomain(const eckit::Parametrisation& p);

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const { return false; }

    static std::string static_type() {return "empty";}
    virtual std::string type() const { return static_type(); }

    virtual bool empty() const { return true; }
    virtual bool global() const { return false; }

    virtual eckit::Properties spec() const;

    virtual void print(std::ostream&) const;

    virtual std::string units() const; // Not implemented

};


}  // namespace domain
}  // namespace atlas
