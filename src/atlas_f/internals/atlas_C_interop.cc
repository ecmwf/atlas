#include <stdlib.h>
#include "eckit/memory/Owned.h"

extern "C"
{
  void atlas_free( void* ptr[] )
  {
    delete[] ptr;
    ptr=NULL;
  }

  int atlas__Owned__owners(const eckit::Owned* owned)
  {
    return owned->owners();
  }

  void atlas__Owned__attach(const eckit::Owned* owned)
  {
    return owned->attach();
  }

  void atlas__Owned__detach(const eckit::Owned* owned)
  {
    return owned->detach();
  }

  int atlas__compare_cptr_equal( void* p1, void* p2 )
  {
    return (p1 == p2);
  }

}



