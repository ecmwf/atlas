
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
atlas::array::ArrayView<typename std::conditional<std::is_const<Value>::value, const byte, byte>::type,Rank+1>
byteView (const atlas::array::ArrayView<Value,Rank> & view)
{
  using B = typename std::conditional<std::is_const<Value>::value, const byte, byte>::type;
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

  B * data = (B *) (view.data ());

  return atlas::array::ArrayView<B,rank> (data, shape, strides);
}

// Add an extra outer dummy dimension

template <typename Value, int Rank>
atlas::array::ArrayView<Value,Rank+1>
addDummyDimension (const atlas::array::ArrayView<Value,Rank> & view)
{
  constexpr int rank = Rank+1;
  atlas::array::ArrayShape shape = {1};
  atlas::array::ArrayStrides strides = {view.size ()};

  for (int i = 0; i < Rank; i++)
    {
      shape.push_back (view.shape (i));
      strides.push_back (view.stride (i));
    }

  using nonConstValue = typename std::remove_const<Value>::type;

  Value * data = (nonConstValue *) (view.data ());

  return atlas::array::ArrayView<Value,rank> (data, shape, strides);
}   



