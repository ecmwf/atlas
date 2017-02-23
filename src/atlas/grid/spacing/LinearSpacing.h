#pragma once

#include <array>
#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

/// @brief Linear spacing in interval
///
/// There are N equally spaced points in the closed interval [start, stop] or the 
/// half-open interval [start, stop) (depending on whether endpoint is True or False)
///
/// Using the constructor LinearSpacing( start, end, N, endpoint ) we can create
///
///    LinearSpacing( 2, 3, 5, true  )                    -->  { 2.0 , 2.25 , 2.5 , 2.75 , 3.0 }
///    LinearSpacing( 2, 3, 5, false )                    -->  { 2.0 , 2.2  , 2.4 , 2.6  , 2.8 }
///
/// Configuration parameters can be passed as well with following keys:
///
///    {"start":2 , "end":3, "N":5, "endpoint":true }     --> { 2.0 , 2.25 , 2.5 , 2.75 , 3.0 }
///    {"start":2 , "end":3, "N":5, "endpoint":false}     --> { 2.0 , 2.2  , 2.4 , 2.6  , 2.8 }
///
/// Instead of the "end" key, you can provide the "length" key, to achieve the same results:
///
///    {"start":2 , "length":1, "N":5, "endpoint":true }  --> { 2.0 , 2.25 , 2.5 , 2.75 , 3.0 }
///    {"start":2 , "length":1, "N":5, "endpoint":false}  --> { 2.0 , 2.2  , 2.4 , 2.6  , 2.8 }


class LinearSpacing: public Spacing {

public:

    /// constructor
    LinearSpacing( const eckit::Parametrisation& p );
    LinearSpacing( double start, double end, long N, bool endpoint=true );

    // LinearSpacing( double centre, double step, long N, bool endpoint=true );

    // class name
    static std::string className() { return "atlas.LinearSpacing"; }
    static std::string spacing_type_str() {return "linear";}

    double step() const;
    
    bool endpoint() const;

public:

    struct Params {
      double start;
      double end;
      long N;
      double length;
      bool endpoint;
      double step;
      Params(const eckit::Parametrisation& p);
    };

protected:

    // points are equally spaced between xmin and xmax
    // Depending on value of endpoint, the spacing will be different
    void setup(double start, double end, long N, bool endpoint);

};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
