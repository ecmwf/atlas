#include <stdlib.h>

extern "C"
{
  void atlas_free(void* ptr[])
  {
    free(*ptr); *ptr=NULL;
  }
}

