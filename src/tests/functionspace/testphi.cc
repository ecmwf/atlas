#include "atlas/field.h"
#include "atlas/array.h"
#include <iostream>

template <typename Value, int Rank>
atlas::array::ArrayView<Value,Rank-1>
dropDimension (const atlas::array::ArrayView<Value,Rank> & view, int dim, atlas::idx_t idx)
{
  constexpr int rank = Rank-1;
  atlas::idx_t _shape[rank], _strides[rank];
  atlas::array::ArrayShape shape;
  atlas::array::ArrayStrides strides;

  if (dim < 0)
    dim = dim + Rank;

  for (int i = 0; i < dim; i++)
    {
      shape.push_back (view.shape (i));
      strides.push_back (view.stride (i));
    }
  atlas::idx_t stride_dim = view.stride (dim);
  for (int i = dim + 1; i < Rank; i++)
    {
      shape.push_back (view.shape (i));
      strides.push_back (view.stride (i));
    }


  using nonConstValue = typename std::remove_const<Value>::type;

  Value * data = (nonConstValue *) (view.data ());

  return atlas::array::ArrayView<Value,rank> (data + stride_dim * idx, shape, strides);
}   

using byte = unsigned char;

#ifdef UNDEF
template <typename Value, int Rank>
atlas::array::ArrayView<byte,Rank+1>
byteView (const atlas::array::ArrayView<Value,Rank> & view)
{
  
}
#endif



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
  atlas::Field f ("f", atlas::array::DataType::kind<double> (), {5, 7});

//auto w = atlas::array::make_view<long,2> (f, atlas::util::Config ("check_type", false));
  auto w = atlas::array::make_view_nocheck<long,2> (f);
  auto v = atlas::array::make_view<double,2> (f);

  prss ("v", v);

  for (int i = 0; i < v.shape ()[0]; i++)
  for (int j = 0; j < v.shape ()[1]; j++)
    v (i, j) = i * 100 + j;


  printf (" &v(0,0) = 0x%llx\n", &v(0,0));
  printf (" &v(1,0) = 0x%llx\n", &v(1,0));
  printf (" &v(0,1) = 0x%llx\n", &v(0,1));

  
  printf (" %8s |", "");
  for (int j = 0; j < v.shape ()[1]; j++)
    printf (" %12d", j);
  printf ("\n");

  for (int i = 0; i < v.shape ()[0]; i++)
    {
      printf (" %8d |", i);
      for (int j = 0; j < v.shape ()[1]; j++)
        printf (" %12.0f", v (i, j));
      printf ("\n");
    }

  
  auto v1 = dropDimension (v, -1, 0);

  prss ("v1", v1);

  std::cout << " v1.rank () = " << v1.rank () << std::endl;

  for (int i = 0; i < v1.shape (0); i++)
    printf (" %8d > %12.2f\n", i, v1 (i));
  

  auto v2 = dropDimension (v, -1, 1);

  std::cout << " v2.rank () = " << v2.rank () << std::endl;

  for (int i = 0; i < v2.shape (0); i++)
    printf (" %8d > %12.2f\n", i, v2 (i));
  
  auto v3 = dropDimension (v, 0, 2);

  std::cout << " v3.rank () = " << v3.rank () << std::endl;

  prss ("v3", v3);

  for (int i = 0; i < v3.shape (0); i++)
    printf (" %8d > %12.2f\n", i, v3 (i));
  

  return 0;

}

