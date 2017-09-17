/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an size_tergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <stddef.h>
#include <vector>

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

class ArrayLayout : public std::vector<size_t> {
private:
   using Base = std::vector<size_t>;
public:
    using Base::Base; // inherit constructors from std::vector<size_t>
};

inline ArrayLayout make_layout(size_t size1) { return ArrayLayout(1,size1); }
inline ArrayLayout make_layout(size_t size1, size_t size2) { ArrayLayout v(2); v[0]=size1; v[1]=size2; return v; }
inline ArrayLayout make_layout(size_t size1, size_t size2, size_t size3) { ArrayLayout v(3); v[0]=size1; v[1]=size2; v[2]=size3; return v; }
inline ArrayLayout make_layout(size_t size1, size_t size2, size_t size3, size_t size4) { ArrayLayout v(4); v[0]=size1; v[1]=size2; v[2]=size3; v[3]=size4; return v; }

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
