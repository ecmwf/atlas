#include "ioFieldDesc.h"
#include "arrayViewHelpers.h"


namespace
{

// The following templates create a list of byte views of rank 2 from a view of any rank, any type
//

template <int Rank>
void listOf1DByteView (atlas::array::ArrayView<byte,Rank> & view,
                       const std::vector<atlas::idx_t> & _ind,
                       const atlas::Field & f, size_t ldim, atlas::idx_t gdim,
                       std::vector<ioFieldDesc> & list)
{
  static_assert (Rank > 2, "listOf1DByteView should be called with views having a Rank > 2");

  std::vector<atlas::idx_t> ind = _ind;
  ind.push_back (0);

  if (gdim == 0)
    {
      for (int i = 0; i < view.shape (1); i++)
        {
          auto v = dropDimension (view, 1, i);
          ind.back () = i;
          listOf1DByteView (v, ind, f, ldim, 0, list);
        }
    }
  else
    {
      for (int i = 0; i < view.shape (0); i++)
        {
          auto v = dropDimension (view, 0, i);
          ind.back () = i;
          listOf1DByteView (v, ind, f, ldim, gdim-1, list);
        }
    }
}

template <>
void listOf1DByteView (atlas::array::ArrayView<byte,2> & view,
                       const std::vector<atlas::idx_t> & ind,
                       const atlas::Field & f, size_t ldim, atlas::idx_t gdim,
                       std::vector<ioFieldDesc> & list)
{
  auto v = addDummyDimension (view);
  list.push_back (ioFieldDesc (v, ind, f, ldim));
}

template <int Rank, typename Value>
void createListOf1DByteView (atlas::array::ArrayView<Value,Rank> & view, 
                             const atlas::Field & f, size_t ldim, atlas::idx_t gdim,
                             std::vector<ioFieldDesc> & list)
{
  auto v = byteView (view);
  listOf1DByteView (v, std::vector<atlas::idx_t> (), f, ldim, gdim, list);
}

};


void createIoFieldDescriptors 
  (atlas::Field & f, std::vector<ioFieldDesc> & list, size_t ldim, atlas::idx_t gdim)
{
  int rank = f.rank ();
  auto type = f.datatype ();

  if (gdim < 0)
    gdim = gdim + rank;

#define HANDLE_TYPE_RANK(__type,__rank) \
   if (rank == __rank)                                                   \
    {                                                                    \
      auto v = atlas::array::make_view<__type,__rank> (f);               \
      createListOf1DByteView (v, f, ldim, gdim, list);                   \
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

void createIoFieldDescriptors 
  (atlas::FieldSet & s, std::vector<ioFieldDesc> & list, size_t ldim, atlas::idx_t gdim)
{
  for (auto & f : s)
    createIoFieldDescriptors (f, list, ldim, gdim);
}



