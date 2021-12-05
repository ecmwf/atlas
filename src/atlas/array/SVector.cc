/*
* (C) Copyright 2013 ECMWF.
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation nor
* does it submit to any jurisdiction.
*/


#include "SVector.h"

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

extern "C" {
void allocate_managedmem_double(double** a, int N) {
    allocate_managedmem(*a, N);
}
void allocate_managedmem_float(double** a, int N) {
    allocate_managedmem(*a, N);
}
void allocate_managedmem_int(int** a, int N) {
    allocate_managedmem(*a, N);
}
void allocate_managedmem_long(long** a, int N) {
    allocate_managedmem(*a, N);
}
}

//------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
