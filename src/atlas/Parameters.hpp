// (C) Copyright 1996-2014 ECMWF.

#ifndef Parameters_hpp
#define Parameters_hpp

namespace atlas {

enum { XX=0, YY=1, ZZ=2 };

struct Entity
{
    enum { NODES=0, FACES=1, ELEMS=2 };
};

struct ElementRef
{
  ElementRef() {}
  ElementRef(int func_space_idx, int elem_idx) :
    f(func_space_idx),e(elem_idx) {}
  int f;
  int e;
};

#if 1

// So that we have the same local id in fortran by natural indexing
#define C_IDX(index) index-1
#define F_IDX(index) index+1

#else

#define C_IDX(index) index
#define F_IDX(index) index

#endif

} // namespace atlas

#endif // Parameters_hpp
