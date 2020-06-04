#include "ioFieldDesc.h"
#include "arrayViewHelpers.h"

#include "atlas/field.h"
#include "atlas/array.h"
#include <iostream>


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

  printf (" shape = ");
  for (int i = 0; i < v.rank (); i++)
    printf (" %8d", v.shape (i));
  printf ("\n");

  for (int jblk = 0; jblk < ngpblks; jblk++)
  for (int jlev = 0, jloc = 0; jlev < nflevg; jlev++)
  for (int jlon = 0; jlon < nproma; jlon++)
    v (jblk, jlev, jlon) = 1000 * jlev + jlon + jblk * nproma;
  
  std::vector<ioFieldDesc> df;
  createIoFieldDescriptorsBlocked (f, df, 0, 2, ngptot); // NGPBLKS dimension, NPROMA dimension, NGPTOT

  for (int i = 0; i < df.size (); i++)
    {
      auto & d = df[i];
      auto & v = d.view ();

      printf (" %8d | %8d | %8d | %8d | %8d\n", i, v.rank (), v.shape (0), v.shape (1), v.shape (2));

      if (i == 3)
        {
          for (int jblk = 0; jblk < v.shape (0); jblk++)
          for (int jlon = 0; jlon < v.shape (1); jlon++)
            {
              long x = 0;
              for (int j = 0; j < 8; j++)
                {
                  byte * b = (byte *)&x + j;
                  *b = v (jblk, jlon, j);
                }
              printf (" %8ld", x);
            }
          printf ("\n");
        }



    }



  return 0;

}

