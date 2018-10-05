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

#include <stddef.h>
#include <vector>
#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

class ArrayAlignment {
public:
    ArrayAlignment() : alignment_( 1 ) {}
    ArrayAlignment( int alignment ) : alignment_( alignment ) {}
    operator int() const { return alignment_; }

private:
    int alignment_;
};

class ArrayShape : public std::vector<idx_t> {
private:
    using Base = std::vector<idx_t>;

public:
    ArrayShape() {}
    ArrayShape( Base&& base ) : Base( std::forward<Base>( base ) ) {}
    ArrayShape( std::initializer_list<idx_t> list ) : Base( list ) {}
};

inline ArrayShape make_shape( std::initializer_list<idx_t> sizes ) {
    return ArrayShape( sizes );
}
inline ArrayShape make_shape( idx_t size1 ) {
    return ArrayShape{size1};
}
inline ArrayShape make_shape( idx_t size1, idx_t size2 ) {
    return ArrayShape{size1, size2};
}
inline ArrayShape make_shape( idx_t size1, idx_t size2, idx_t size3 ) {
    return ArrayShape{size1, size2, size3};
}
inline ArrayShape make_shape( idx_t size1, idx_t size2, idx_t size3, idx_t size4 ) {
    return ArrayShape{size1, size2, size3, size4};
}
inline ArrayShape make_shape( idx_t size1, idx_t size2, idx_t size3, idx_t size4, idx_t size5 ) {
    return ArrayShape{size1, size2, size3, size4, size5};
}

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
