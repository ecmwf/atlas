#include <stdlib.h>
extern "C"
{
  void atlas_free( void* ptr[] )
  {
    delete[] ptr;
    ptr=NULL;
  }
}

