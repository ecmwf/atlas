#ifndef Parameters_hpp
#define Parameters_hpp

namespace ecmwf {

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

} // namespace ecmwf

#endif // Parameters_hpp
