#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

#include "FloatCompare.h"

using namespace std;

//------------------------------------------------------------------------------------------------------

namespace eckit {

//------------------------------------------------------------------------------------------------------

// See:
// * http://randomascii.wordpress.com/2012/01/11/tricks-with-the-floating-point-format
// * http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition
//   for the potential portability problems with the union and bit-fields below.

/// @TODO check in BIGENDIAN archs...

union Float_t
{
    Float_t(float num = 0.0f) : f(num) { assert( sizeof(i) == sizeof(f) ); }

    // Portable extraction of components.
    bool Negative() const { return (i >> 31) != 0; }
    int32_t RawMantissa() const { return i & ((1 << 23) - 1); }
    int32_t RawExponent() const { return (i >> 23) & 0xFF; }

    int32_t i;
    float f;
};

union Double_t
{
    Double_t(double num = 0.0) : f(num) { assert( sizeof(i) == sizeof(f) ); }
    // Portable extraction of components.
    bool Negative() const { return (i >> 63) != 0; }
    int64_t RawMantissa() const { return i & ((1LL << 52) - 1); }
    int64_t RawExponent() const { return (i >> 52) & 0x7FF; }

    int64_t i;
    double f;
};

template < class T > struct FPCompare;

template <> struct FPCompare<float>  { typedef Float_t  FP_t; };
template <> struct FPCompare<double> { typedef Double_t FP_t; };

//------------------------------------------------------------------------------------------------------

template < typename T >
bool AlmostEqualUlps(T A, T B, int maxUlpsDiff)
{
    typename FPCompare<T>::FP_t uA(A);
    typename FPCompare<T>::FP_t uB(B);

    // Different signs means they do not match.
    if (uA.Negative() != uB.Negative())
    {
        // Check for equality to make sure +0==-0
        if (A == B)
            return true;
        return false;
    }

    // Find the difference in ULPs.
    int ulpsDiff = abs(int(uA.i - uB.i));
    if (ulpsDiff <= maxUlpsDiff)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------

template < typename T >
bool AlmostEqualUlpsAndAbs(T A, T B, T maxDiff, int maxUlpsDiff)
{
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    float absDiff = fabs(A - B);
    if (absDiff <= maxDiff)
        return true;

    typename FPCompare<T>::FP_t uA(A);
    typename FPCompare<T>::FP_t uB(B);

    // Different signs means they do not match.
    if (uA.Negative() != uB.Negative())
        return false;

    // Find the difference in ULPs.
    int ulpsDiff = abs(int(uA.i - uB.i));
    if (ulpsDiff <= maxUlpsDiff)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------

template < typename T >
bool AlmostEqualRelativeAndAbs(T A, T B, T maxDiff, T maxRelDiff)
{
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    T diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;

    A = fabs(A);
    B = fabs(B);
    T largest = (B > A) ? B : A;

    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}

//-----------------------------------------------------------------------------

bool FloatCompare::is_equal(double a, double b)
{
//    return AlmostEqualUlps(a,b,10);
    return AlmostEqualUlpsAndAbs(a,b,std::numeric_limits<double>::epsilon(),10);
}

//-----------------------------------------------------------------------------

} // namespace eckit
