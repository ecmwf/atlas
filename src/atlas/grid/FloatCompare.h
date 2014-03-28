#ifndef eckit_FloatCompare_h
#define eckit_FloatCompare_h

//-----------------------------------------------------------------------------

namespace eckit {

//-----------------------------------------------------------------------------

/// @todo move this class into eckit

class FloatCompare {
public:

    static bool is_equal( float  a, float  b );
    static bool is_equal( double a, double b );

};

//---------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif

