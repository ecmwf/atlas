/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Base64.h"

#include <arpa/inet.h>
#include <array>

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

namespace {

class Base64Tables {
public:
    std::array<unsigned char, 256> decode;
    std::string encode{
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/"};
    static const Base64Tables& instance() {
        static Base64Tables instance;
        return instance;
    }

private:
    Base64Tables() {
        std::array<unsigned char, 256> tmp;

        const char* p = encode.data();

        size_t i = 0;
        while ( *p ) {
            size_t j  = *p;
            tmp[i]    = *p;
            decode[j] = i;
            p++;
            i++;
        }
    }
};

}  // namespace

//---------------------------------------------------------------------------------------------------------------------

std::string Base64::encode( const void* data, size_t len ) {
    const auto& table        = Base64Tables::instance().encode;
    const unsigned char* src = reinterpret_cast<const unsigned char*>( data );
    unsigned char *out, *pos;
    const unsigned char *end, *in;

    size_t olen;

    olen = 4 * ( ( len + 2 ) / 3 ); /* 3-byte blocks to 4-byte */

    if ( olen < len )
        return std::string(); /* integer overflow */

    std::string outStr;
    outStr.resize( olen );
    out = (unsigned char*)&outStr[0];

    end = src + len;
    in  = src;
    pos = out;
    while ( end - in >= 3 ) {
        *pos++ = table[in[0] >> 2];
        *pos++ = table[( ( in[0] & 0x03 ) << 4 ) | ( in[1] >> 4 )];
        *pos++ = table[( ( in[1] & 0x0f ) << 2 ) | ( in[2] >> 6 )];
        *pos++ = table[in[2] & 0x3f];
        in += 3;
    }

    if ( end - in ) {
        *pos++ = table[in[0] >> 2];
        if ( end - in == 1 ) {
            *pos++ = table[( in[0] & 0x03 ) << 4];
            *pos++ = '=';
        }
        else {
            *pos++ = table[( ( in[0] & 0x03 ) << 4 ) | ( in[1] >> 4 )];
            *pos++ = table[( in[1] & 0x0f ) << 2];
        }
        *pos++ = '=';
    }
    return outStr;
}

//---------------------------------------------------------------------------------------------------------------------

std::string Base64::decode( const void* data, size_t len ) {
    const auto& table = Base64Tables::instance().decode;

    unsigned char* p = (unsigned char*)data;
    int pad          = len > 0 && ( len % 4 || p[len - 1] == '=' );
    const size_t L   = ( ( len + 3 ) / 4 - pad ) * 4;
    std::string str( L / 4 * 3 + pad, '\0' );

    for ( size_t i = 0, j = 0; i < L; i += 4 ) {
        int n    = table[p[i]] << 18 | table[p[i + 1]] << 12 | table[p[i + 2]] << 6 | table[p[i + 3]];
        str[j++] = n >> 16;
        str[j++] = n >> 8 & 0xFF;
        str[j++] = n & 0xFF;
    }
    if ( pad ) {
        int n               = table[p[L]] << 18 | table[p[L + 1]] << 12;
        str[str.size() - 1] = n >> 16;

        if ( len > L + 2 && p[L + 2] != '=' ) {
            n |= table[p[L + 2]] << 6;
            str.push_back( n >> 8 & 0xFF );
        }
    }
    return str;
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
