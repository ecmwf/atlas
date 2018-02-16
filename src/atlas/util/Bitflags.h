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

#include <string>

namespace atlas {
namespace util {

class Bitflags {
public:
    static void reset( int& flags, int bit = 0 ) { flags = bit; }

    static void set( int& flags, int bit ) { flags |= bit; }

    static void unset( int& flags, int bit ) { flags &= ( ~bit ); }

    static void toggle( int& flags, int bit ) { flags ^= bit; }

    static bool check( int flags, int bits ) { return ( flags & bits ) == bits; }

    static bool check_all( int flags, int bits ) { return ( flags & bits ) == bits; }

    static bool check_any( int flags, int bits ) { return flags & bits; }

    static std::string bitstr( int flags ) {
        char str[9] = {0};
        int i;
        for ( i = 7; i >= 0; i-- ) {
            str[i] = ( flags & 1 ) ? '1' : '0';
            flags >>= 1;
        }
        return std::string( str, 9 );
    }
};

}  // namespace util
}  // namespace atlas
