#include "atlas/field.h"
#include "atlas/array.h"
#include <iostream>



int main (int argc, char * argv[])
{
  atlas::Field f ("f", atlas::array::DataType::kind<double> (), {5, 7});

//auto w = atlas::array::make_view<long,2> (f, atlas::util::Config ("check_type", false));
  auto w = atlas::array::make_view_nocheck<long,2> (f);
  auto v = atlas::array::make_view<double,2> (f);

  for (int i = 0; i < v.shape ()[0]; i++)
  for (int j = 0; j < v.shape ()[1]; j++)
    v (i, j) = i * 100 + j;
  

  auto v1 = v.drop (0);

  std::cout << " v1.rank () = " << v1.rank () << std::endl;

  for (int i = 0; i < v1.shape (0); i++)
    printf (" %8d > %12.2f\n", i, v1 (i));
  

  auto v2 = v.drop (1);

  std::cout << " v2.rank () = " << v2.rank () << std::endl;

  for (int i = 0; i < v2.shape (0); i++)
    printf (" %8d > %12.2f\n", i, v2 (i));
  

  return 0;

}

