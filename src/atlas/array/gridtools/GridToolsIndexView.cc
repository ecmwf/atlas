/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <stdexcept>
#include "atlas/field/Field.h"
#include "atlas/array/Array.h"
#include "atlas/array/IndexView.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void IndexView<Value,Rank>::dump(std::ostream& os) const {
  os << "size: " << size() << " , values: ";
  os << "[ ";
  for( size_t j=0; j<size(); ++ j )
    os << (*this)(j) << " ";
  os << "]" << std::endl;
}

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas

//------------------------------------------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {

template class IndexView<int,1>;

}
}
