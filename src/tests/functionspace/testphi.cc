#include "atlas/field.h"
#include "atlas/array.h"
#include <iostream>


// Drop a dimension of a view (set it to a fixed value); we get a view of rank=Rank-1

template <typename Value, int Rank>
atlas::array::ArrayView<Value,Rank-1>
dropDimension (const atlas::array::ArrayView<Value,Rank> & view, int dim, atlas::idx_t idx)
{
  constexpr int rank = Rank-1;
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


// Convert a view to a view of bytes (add an extra inner dimension)

using byte = unsigned char;

template <typename Value, int Rank>
atlas::array::ArrayView<byte,Rank+1>
byteView (const atlas::array::ArrayView<Value,Rank> & view)
{
  constexpr int rank = Rank+1;
  atlas::array::ArrayShape shape;
  atlas::array::ArrayStrides strides;

  size_t dlen = atlas::array::DataType::create<Value> ().size ();

  for (int i = 0; i < Rank; i++)
    {
      shape.push_back (view.shape (i));
      strides.push_back (view.stride (i) * dlen);
    }

  shape.push_back (dlen);
  strides.push_back (1);

  byte * data = (byte *) (view.data ());

  return atlas::array::ArrayView<byte,rank> (data, shape, strides);
}


// The following templates create a list of byte views of rank 2 from a view of any rank, any type
//

class ioFieldDesc
{
public:
  ioFieldDesc (atlas::array::ArrayView<byte,2> & _v, 
               const std::vector<atlas::idx_t> & _ind,
               const atlas::Field & _f) : v (_v), ind (_ind), f (_f) {}
  const std::vector<atlas::idx_t> & indices () const
  {
    return ind;
  }
  const atlas::Field & field () const
  {
    return f;
  }
private:
  atlas::array::ArrayView<byte,2> v;
  const std::vector<atlas::idx_t> ind;
  const atlas::Field & f;
};

template <int Rank>
void listOf1DByteView (atlas::array::ArrayView<byte,Rank> & view,
                       const std::vector<atlas::idx_t> & _ind,
                       const atlas::Field & f,
                       std::vector<ioFieldDesc> & list)
{
  static_assert (Rank > 2, "listOf1DByteView should be called with views having a Rank > 2");

  std::vector<atlas::idx_t> ind = _ind;
  ind.push_back (0);

  for (int i = 0; i < view.shape (0); i++)
    {
      auto v = dropDimension (view, 0, i);
      ind.back () = i;
      listOf1DByteView (v, ind, f, list);
    }
}

template <>
void listOf1DByteView (atlas::array::ArrayView<byte,2> & view,
                       const std::vector<atlas::idx_t> & ind,
                       const atlas::Field & f,
                       std::vector<ioFieldDesc> & list)
{
  list.push_back (ioFieldDesc (view, ind, f));
}

template <int Rank, typename Value>
void createListOf1DByteView (atlas::array::ArrayView<Value,Rank> & view, 
                             const atlas::Field & f,
                             std::vector<ioFieldDesc> & list)
{
  auto v = byteView (view);
  listOf1DByteView (v, std::vector<atlas::idx_t> (), f, list);
}

void createIoFieldDescriptors (atlas::Field & f, std::vector<ioFieldDesc> & list)
{
  int rank = f.rank ();
  auto type = f.datatype ();

#define HANDLE_TYPE_RANK(__type,__rank) \
   if (rank == __rank)                                                   \
    {                                                                    \
      auto v = atlas::array::make_view<__type,__rank> (f);               \
      createListOf1DByteView (v, f, list);                               \
      goto done;                                                         \
    }

#define HANDLE_TYPE(__type) \
  if (type.kind () == atlas::array::DataType::create<__type> ().kind ()) \
    {                                                                                           \
      HANDLE_TYPE_RANK (__type, 1); HANDLE_TYPE_RANK (__type, 2); HANDLE_TYPE_RANK (__type, 3); \
      HANDLE_TYPE_RANK (__type, 4); HANDLE_TYPE_RANK (__type, 5); HANDLE_TYPE_RANK (__type, 6); \
      HANDLE_TYPE_RANK (__type, 7); HANDLE_TYPE_RANK (__type, 8); HANDLE_TYPE_RANK (__type, 9); \
    }

  HANDLE_TYPE (long);
  HANDLE_TYPE (double);
  HANDLE_TYPE (int);
  HANDLE_TYPE (float);

done:

  return;
}


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

