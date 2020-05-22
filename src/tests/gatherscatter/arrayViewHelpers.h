
#pragma once


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

