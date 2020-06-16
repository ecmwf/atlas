#include "ioFieldDesc.h"
#include "arrayViewHelpers.h"

#include "atlas/field.h"
#include "atlas/array.h"
#include <iostream>
#include <stdio.h>


template <typename V>
void prss (const std::string & t, const V & v)
{
  printf ("--------- %s ---------\n", t.data ());
  for (int i = 0; i < v.rank (); i++)
    printf (" %8d > %8d, %8d\n", i, v.shape (i), v.stride (i));
  printf ("\n");
}

int main (int argc, char * argv[])
{
  const int ngpblks = 10, nflevg = 20, nproma = 16;

  const int ngptot = ngpblks * nproma;

  atlas::Field f ("f", atlas::array::DataType::kind<long> (), {ngpblks, nflevg, nproma});

  auto v = atlas::array::make_view<long,3> (f);

  viewLoop (v, [] (long & x, int i, int j, int k) { x = 100000 * i + 100 * j + k; });

  for (int i = 0; i < v.shape (0); i++)
  for (int j = 0; j < v.shape (1); j++)
  for (int k = 0; k < v.shape (2); k++)
    printf (" %10.10d\n", v (i, j, k));

  return 0;

}

