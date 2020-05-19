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
  atlas::Field f ("f", atlas::array::DataType::kind<long> (), {5, 7});

  auto v = atlas::array::make_view<long,2> (f);

  prss ("v", v);

  for (int i = 0; i < v.shape ()[0]; i++)
  for (int j = 0; j < v.shape ()[1]; j++)
    v (i, j) = i * 100 + j;


  printf (" &v(0,0) = 0x%llx\n", &v(0,0));
  printf (" &v(1,0) = 0x%llx\n", &v(1,0));
  printf (" &v(0,1) = 0x%llx\n", &v(0,1));

  auto c = byteView (v);
  
  std::vector<ioFieldDesc> list;

  createIoFieldDescriptors (f, list);

  std::cout << " list.size () = " << list.size () << std::endl;

  prss ("c", c);

  printf (" sizeof (c (0,0)) = %d\n", sizeof (c (0,0)));

  for (int i = 0; i < 4; i++)
  printf (" &c(0,0,%d) = 0x%llx\n", i, &c(0,0,i));

  for (int i = 0; i < 4; i++)
  printf (" &c(0,%d,0) = 0x%llx\n", i, &c(0,i,0));

  for (int i = 0; i < 4; i++)
  printf (" &c(%d,0,0) = 0x%llx\n", i, &c(i,0,0));


if(1)
  for (int i = 0; i < 3; i++)
    {
      c (i, i, 0) = 0xff;
      c (i, i, 1) = 0xff;
      c (i, i, 2) = 0xff;
      c (i, i, 3) = 0xff;
      c (i, i, 4) = 0xff;
      c (i, i, 5) = 0xff;
      c (i, i, 6) = 0xff;
      c (i, i, 7) = 0xff;
    }
  
  printf (" %8s |", "");
  for (int j = 0; j < v.shape ()[1]; j++)
    printf (" %12d", j);
  printf ("\n");

  for (int i = 0; i < v.shape ()[0]; i++)
    {
      printf (" %8d |", i);
      for (int j = 0; j < v.shape ()[1]; j++)
        printf (" %12d", v (i, j));
      printf ("\n");
    }

  
  auto v1 = dropDimension (v, -1, 0);

  prss ("v1", v1);

  std::cout << " v1.rank () = " << v1.rank () << std::endl;

  for (int i = 0; i < v1.shape (0); i++)
    printf (" %8d > %12d\n", i, v1 (i));
  

  auto v2 = dropDimension (v, -1, 1);

  std::cout << " v2.rank () = " << v2.rank () << std::endl;

  for (int i = 0; i < v2.shape (0); i++)
    printf (" %8d > %12d\n", i, v2 (i));
  
  auto v3 = dropDimension (v, 0, 2);

  std::cout << " v3.rank () = " << v3.rank () << std::endl;

  prss ("v3", v3);

  for (int i = 0; i < v3.shape (0); i++)
    printf (" %8d > %12d\n", i, v3 (i));
  

  return 0;

}

