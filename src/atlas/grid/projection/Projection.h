#pragma once

#include "eckit/geometry/KPoint.h"
#include "eckit/geometry/Point2.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/value/Properties.h"
#include "eckit/memory/Builder.h"
#include "eckit/memory/Owned.h"

namespace atlas {
namespace grid {
namespace projection {

class Projection : public eckit::Owned {

public:

    typedef const eckit::Parametrisation& ARG1;
    typedef eckit::BuilderT1<Projection> builder_t;

public:

    static Projection* create(); // creates the LonLatProjection
    static Projection* create(const eckit::Parametrisation& p);

    Projection() {}
    Projection( const Projection& rhs ) {}   // copy constructor
    virtual Projection* clone() const =0;    // clone method acting like virtual copy constructor
    virtual ~Projection() {}                 // destructor should be virtual when using a virtual copy constructor

    static std::string className() {return "atlas.Projection";}
    static std::string projection_type_str() {return "projection";}
    virtual std::string virtual_projection_type_str() const { return "projection"; }

    // purely virtual functions: must be implemented by inheriting classes
    virtual eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2) const =0;
    virtual eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2) const =0;
    virtual bool isRegional()=0;

    virtual eckit::Properties spec() const =0;

};

}  // namespace projection
}  // namespace grid
}  // namespace atlas
