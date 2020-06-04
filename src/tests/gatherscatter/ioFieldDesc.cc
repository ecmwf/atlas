#include "ioFieldDesc.h"
#include "arrayViewHelpers.h"
#include "atlas/runtime/Exception.h"


namespace
{

// The following templates create a list of byte views of rank 2 from a view of any rank, any type


// Non-blocked (1:NGPTOT)

template <int Rank>
void listOf1DByteView (atlas::array::ArrayView<byte,Rank> & view,
                       const std::vector<atlas::idx_t> & _ind,
                       const atlas::Field & f, size_t ngptot, atlas::idx_t gdim,
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
          listOf1DByteView (v, ind, f, ngptot, 0, list);
        }
    }
  else
    {
      for (int i = 0; i < view.shape (0); i++)
        {
          auto v = dropDimension (view, 0, i);
          ind.back () = i;
          listOf1DByteView (v, ind, f, ngptot, gdim-1, list);
        }
    }
}

template <>
void listOf1DByteView (atlas::array::ArrayView<byte,2> & view,
                       const std::vector<atlas::idx_t> & ind,
                       const atlas::Field & f, size_t ngptot, atlas::idx_t gdim,
                       std::vector<ioFieldDesc> & list)
{
  auto v = addDummyDimension (view);
  list.push_back (ioFieldDesc (v, ind, f, ngptot));
}

template <int Rank, typename Value>
void createListOf1DByteView (atlas::array::ArrayView<Value,Rank> & view, 
                             const atlas::Field & f, size_t ngptot, atlas::idx_t gdim,
                             std::vector<ioFieldDesc> & list)
{
  auto v = byteView (view);
  listOf1DByteView (v, std::vector<atlas::idx_t> (), f, ngptot, gdim, list);
}

// Blocked (1:NPROMA,1:NGPBLKS)

template <int Rank>
void listOf1DByteViewBlocked 
  (atlas::array::ArrayView<byte,Rank> & view,
   const std::vector<atlas::idx_t> & _ind,
   const atlas::Field & f, atlas::idx_t bdim, atlas::idx_t gdim, size_t ngptot,
   std::vector<ioFieldDesc> & list)
{
  static_assert (Rank > 3, "listOf1DByteViewBlocked should be called with views having a Rank > 3");

  std::vector<atlas::idx_t> ind = _ind;
  ind.push_back (0);


  // First dimension is block
  if (bdim == 0)
    {
      // Grid dimension is second
      if (gdim == 1)
        {
          // Drop dimensions after grid dimension
          for (int i = 0; i < view.shape (2); i++)
            {
              auto v = dropDimension (view, 2, i);
              ind.back () = i;
              listOf1DByteViewBlocked (v, ind, f, bdim, gdim, ngptot, list);
            }
        }
      else
        {
          // Remove second dimension
          for (int i = 0; i < view.shape (1); i++)
            {
              auto v = dropDimension (view, 1, i);
              ind.back () = i;
              listOf1DByteViewBlocked (v, ind, f, bdim, gdim-1, ngptot, list);
            }
        }
    }
  else
    {
      // Remove another leading dimension
      for (int i = 0; i < view.shape (0); i++)
        {
          auto v = dropDimension (view, 0, i);
          ind.back () = i;
          listOf1DByteViewBlocked (v, ind, f, bdim-1, gdim-1, ngptot, list);
        }
    }

}

template <>
void listOf1DByteViewBlocked 
  (atlas::array::ArrayView<byte,3> & view,
   const std::vector<atlas::idx_t> & ind,
   const atlas::Field & f, atlas::idx_t bdim, atlas::idx_t gdim, size_t ngptot, 
   std::vector<ioFieldDesc> & list)
{
  list.push_back (ioFieldDesc (view, ind, f, ngptot));
}

template <int Rank, typename Value>
void createListOf1DByteViewBlocked (atlas::array::ArrayView<Value,Rank> & view, 
                                    const atlas::Field & f, atlas::idx_t bdim, atlas::idx_t gdim, 
                                    size_t ngptot, std::vector<ioFieldDesc> & list)
{
  auto v = byteView (view);
  listOf1DByteViewBlocked (v, std::vector<atlas::idx_t> (), f, bdim, gdim, ngptot, list);
}

};


void createIoFieldDescriptorsBlocked
  (atlas::Field & f, std::vector<ioFieldDesc> & list, atlas::idx_t bdim, atlas::idx_t gdim, size_t ngptot)
{
  int rank = f.rank ();
  auto type = f.datatype ();
  
  if (gdim < 0)
    gdim = gdim + rank;

#define HANDLE_TYPE_RANK(__type,__rank) \
   if (rank == __rank)                                                   \
    {                                                                    \
      auto v = atlas::array::make_view<__type,__rank> (f);               \
      createListOf1DByteViewBlocked (v, f, bdim, gdim, ngptot, list);    \
      goto done;                                                         \
    }

#define HANDLE_TYPE(__type) \
  if (type.kind () == atlas::array::DataType::create<__type> ().kind ()) \
    {                                                                                           \
                                    HANDLE_TYPE_RANK (__type, 2); HANDLE_TYPE_RANK (__type, 3); \
      HANDLE_TYPE_RANK (__type, 4); HANDLE_TYPE_RANK (__type, 5); HANDLE_TYPE_RANK (__type, 6); \
      HANDLE_TYPE_RANK (__type, 7); HANDLE_TYPE_RANK (__type, 8); HANDLE_TYPE_RANK (__type, 9); \
    }

  HANDLE_TYPE (long);
  HANDLE_TYPE (double);
  HANDLE_TYPE (int);
  HANDLE_TYPE (float);

#undef HANDLE_TYPE_RANK
#undef HANDLE_TYPE

  atlas::throw_NotImplemented ("createIoFieldDescriptorsBlocked type/rank", Here ());

done:

  return;
}


void createIoFieldDescriptors 
  (atlas::Field & f, std::vector<ioFieldDesc> & list, size_t ngptot, atlas::idx_t gdim)
{
  int rank = f.rank ();
  auto type = f.datatype ();

  if (gdim < 0)
    gdim = gdim + rank;

#define HANDLE_TYPE_RANK(__type,__rank) \
   if (rank == __rank)                                                   \
    {                                                                    \
      auto v = atlas::array::make_view<__type,__rank> (f);               \
      createListOf1DByteView (v, f, ngptot, gdim, list);                 \
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

#undef HANDLE_TYPE_RANK
#undef HANDLE_TYPE

  atlas::throw_NotImplemented ("createIoFieldDescriptors type/rank", Here ());

done:

  return;
}

void createIoFieldDescriptors 
  (atlas::FieldSet & s, std::vector<ioFieldDesc> & list, size_t ngptot, atlas::idx_t gdim)
{
  for (auto & f : s)
    createIoFieldDescriptors (f, list, ngptot, gdim);
}



