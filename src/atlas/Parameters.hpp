// (C) Copyright 1996-2014 ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation nor
// does it submit to any jurisdiction.


#ifndef Parameters_hpp
#define Parameters_hpp

namespace atlas {

enum { XX=0, YY=1, ZZ=2 };
enum { NODES=0, FACES=1, ELEMS=2 };

struct ElementRef
{
  ElementRef() {}
  ElementRef(int func_space_idx, int elem_idx) :
    f(func_space_idx),e(elem_idx) {}
  int f;
  int e;
};

// So that we have the same local id in fortran by natural indexing
#define C_IDX(index) index-1
#define F_IDX(index) index+1


} // namespace atlas

#endif // Parameters_hpp
