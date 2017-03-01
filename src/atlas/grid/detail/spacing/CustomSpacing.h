#pragma once

#include <array>
#include "atlas/grid/detail/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

/// @brief Custom spacing in interval
///
/// There are N points, by default in the open interval (90, -90).
///
/// Using the constructor CustomSpacing( N, x[], {min,max} ) we can create
///
///    CustomSpacing( 4, {75,25,-25,-75} )                  -->  { 75 , 25 , -25 , -75 }
///    CustomSpacing( 4, {75,25,-25,-75}, {-90,90} )        -->  { 75 , 25 , -25 , -75 }
///
/// The optional argument {min,max} serves as purpose to indicate that the points
/// lie in the open interval (min,max). If not specified, the default values are taken
/// to be the North and South pole's latitudes.
///
/// Configuration parameters can be passed as well with following keys:
///
///    {"N":4, "values":[75,25,-25,75] }                        -->  { 75 , 25 , -25 , -75 }
///    {"N":4, "values":[75,25,-25,75], "interval":[-90,90] }   -->  { 75 , 25 , -25 , -75 }

class CustomSpacing: public Spacing {
private:
    using Interval = std::array<double,2>;
    static constexpr double North() { return  90.; }
    static constexpr double South() { return -90.; }
public:

    // constructor
    CustomSpacing(const eckit::Parametrisation& p);

    CustomSpacing(long N, const double x[], const Interval& = {North(),South()} );

    // class name
    static std::string static_type() {return "custom";}
    virtual std::string type() const {return static_type();}
    
    virtual eckit::Properties spec() const;

};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
