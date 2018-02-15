#ifndef atlas_GaussianSpacing_H
#define atlas_GaussianSpacing_H

#include "atlas/grid/detail/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

/// @brief Gaussian spacing in interval
///
/// There are N Gaussian spaced points in the open interval (90, -90)
///
/// Using the constructor GaussianSpacing( N ) we can create
///
///    Gaussian( 4 )                       -->  { 59.44... , 19.87... ,
///    -19.87... , -59.44... }
///
/// Configuration parameters can be passed as well with following keys:
///
///    {"N":4 }                            -->  { 59.44... , 19.87... ,
///    -19.87... , -59.44... }
///
/// To reverse the orientation of points to go from negative to positive
/// instead, pass also
/// the start and end keys:
///
///    {"N":4, "start":-90, "end":90 }    -->  { -59.44... , -19.87... ,
///    19.87... , 59.44... }

class GaussianSpacing : public Spacing {
public:
    // constructor
    GaussianSpacing( const eckit::Parametrisation& p );

    GaussianSpacing( long N );

    // class name
    static std::string static_type() { return "gaussian"; }
    virtual std::string type() const { return static_type(); }

    virtual Spec spec() const;
};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

#endif
