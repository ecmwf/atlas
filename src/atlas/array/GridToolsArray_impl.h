#pragma once

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

/// TODO: make_view to forward declarations
template <typename Value, unsigned int NDims, bool ReadOnly = false>
static ArrayView<Value, NDims>
make_view(const Array& array);

template <typename Value, unsigned int RANK, unsigned int Dim>
struct array_initializer_impl {

  static void apply(Array const& orig, Array& array_resized) {
      array_initializer_impl<Value, RANK, Dim>::apply(make_view<Value, RANK>(orig), make_view<Value, RANK>(array_resized));
  }
  template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& array_resized, DimIndex... idxs) {
      for(size_t i=0; i < orig.shape(Dim); ++i)
      {
          array_initializer_impl<Value, RANK, Dim+1>::apply(std::move(orig), std::move(array_resized), idxs..., i);
      }
  }
};

template <typename Value, unsigned int RANK>
struct array_initializer_impl<Value, RANK, RANK> {
    template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& array_resized, DimIndex... idxs) {
      array_resized(idxs...) = orig(idxs...);
  }
};

template <unsigned int RANK>
struct array_initializer {
  template <typename... DimIndex>
  static void apply(Array const& orig, Array& array_resized, DimIndex... idxs) {
    switch (orig.datatype().kind()) {
      case DataType::KIND_REAL64:
        array_initializer_impl<double, RANK, 0>::apply(orig, array_resized, idxs...);
        break;
      case DataType::KIND_REAL32:
        array_initializer_impl<float, RANK, 0>::apply(orig, array_resized, idxs...);
        break;
      case DataType::KIND_INT32:
        array_initializer_impl<int, RANK, 0>::apply(orig, array_resized, idxs...);
        break;
      case DataType::KIND_INT64:
        array_initializer_impl<long, RANK, 0>::apply(orig, array_resized, idxs...);
        break;
        ;
      default: {
        std::stringstream err;
        err << "data kind " << orig.datatype().kind() << " not recognised.";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }
};

template<typename Value, unsigned int RANK, unsigned int Dim, unsigned int PartDim>
struct array_initializer_partitioned_val_impl {
  static void apply(Array const& orig, Array& dest, unsigned int pos) {
      auto view = make_view<Value, RANK>(orig);
      array_initializer_partitioned_val_impl<Value, RANK, Dim, PartDim>::apply(make_view<Value, RANK>(orig), make_view<Value, RANK>(dest), pos);
  }
  template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& dest, unsigned int pos, DimIndex... idxs) {
      for(size_t i=0; i < orig.shape(Dim); ++i)
      {
          unsigned int offset = i;
          if(Dim == PartDim && i >= pos) offset += pos;

          std::pair<int,int> pair_idx{i,offset};
          array_initializer_partitioned_val_impl<Value, RANK, Dim+1, PartDim>::apply(std::move(orig), std::move(dest), pos, idxs..., pair_idx);
      }
  }

};

template <typename Value, unsigned int RANK, unsigned int PartDim>
struct array_initializer_partitioned_val_impl<Value, RANK, RANK, PartDim> {
    template <typename ... DimIndex>
  static void apply(ArrayView<Value, RANK> const&& orig, ArrayView<Value, RANK>&& dest, unsigned int pos, DimIndex... idxs) {
      dest(std::get<1>(idxs)...) = orig(std::get<0>(idxs)...);
  }
};


template<unsigned int RANK, unsigned int PartDim>
struct array_initializer_partitioned_impl {
  static void apply(Array const& orig, Array& dest, unsigned int pos) {
      switch (orig.datatype().kind()) {
        case DataType::KIND_REAL64:
          array_initializer_partitioned_val_impl<double, RANK, 0, PartDim>::apply(orig, dest, pos);
          break;
        case DataType::KIND_REAL32:
          array_initializer_partitioned_val_impl<float, RANK, 0, PartDim>::apply(orig, dest,pos );
          break;
        case DataType::KIND_INT32:
          array_initializer_partitioned_val_impl<int, RANK, 0, PartDim>::apply(orig, dest, pos);
          break;
        case DataType::KIND_INT64:
          array_initializer_partitioned_val_impl<long, RANK, 0, PartDim>::apply(orig, dest, pos);
          break;
          ;
        default: {
          std::stringstream err;
          err << "data kind " << orig.datatype().kind() << " not recognised.";
          throw eckit::BadParameter(err.str(), Here());
        }
      }
  }
};

template<unsigned int PartDim>
struct array_initializer_partitioned {
  static void apply(Array const& orig, Array& dest, unsigned int pos) {
    switch (orig.rank()) {
      case 1:
        array_initializer_partitioned_impl<1, PartDim>::apply(orig, dest, pos);
        break;
      case 2:
        array_initializer_partitioned_impl<2, PartDim>::apply(orig, dest, pos);
        break;
      case 3:
        array_initializer_partitioned_impl<3, PartDim>::apply(orig, dest, pos);
        break;
      case 4:
        array_initializer_partitioned_impl<4, PartDim>::apply(orig, dest, pos);
        break;
      case 5:
        array_initializer_partitioned_impl<5, PartDim>::apply(orig, dest, pos);
        break;
      case 6:
        array_initializer_partitioned_impl<6, PartDim>::apply(orig, dest, pos);
        break;
      case 7:
        array_initializer_partitioned_impl<7, PartDim>::apply(orig, dest, pos);
        break;
      case 8:
        array_initializer_partitioned_impl<8, PartDim>::apply(orig, dest, pos);
        break;
      case 9:
        array_initializer_partitioned_impl<9, PartDim>::apply(orig, dest, pos);
        break;
      default: {
        std::stringstream err;
        err << "too high rank";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }
};

inline void Array::insert(size_t idx1, size_t size1)
{
    if(!data_store_->is_on_host()) {
        data_store_->clone_from_device();
    }

    ArrayShape nshape = shape();
    nshape[0] += size1;
    Array* array_resized = Array::create(datatype(), nshape);

    array_initializer_partitioned<0>::apply( *this, *array_resized, idx1);
    data_store_.swap(array_resized->data_store());

    spec_ = array_resized->spec();
//TODO when deleting this if seg fault
//    delete array_resized;

}

}
}
