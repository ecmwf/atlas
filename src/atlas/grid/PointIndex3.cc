#include "PointIndex3.h"
#include "FloatCompare.h"

//------------------------------------------------------------------------------------------------------

using namespace atlas;

namespace eckit {

//------------------------------------------------------------------------------------------------------

bool points_equal(const KPoint3 &a, const KPoint3 &b)
{
    return FloatCompare::is_equal( eckit::KPoint3::distance2(a,b), 0.0 );
}

//------------------------------------------------------------------------------------------------------

} // namespace eckit

