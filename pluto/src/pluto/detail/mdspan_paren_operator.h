/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// This file is to be included from pluto/mdspan.h
// It defines a derived class of std::mdspan with paren operator to access indices.

namespace detail{
  template<class IndexType, class ... Arguments>
  inline constexpr bool valid_indices() {
    return
        ((std::is_convertible<Arguments, IndexType>::value) && ...) &&
        ((std::is_nothrow_constructible<IndexType, Arguments>::value) && ...);
}
}

template<typename ElementType, typename Extents, typename LayoutPolicy = layout_right, typename AccessorPolicy = default_accessor<ElementType>>
class mdspan : public ::STD_MDSPAN_NAMESPACE::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy> {
    using base_type_ = ::STD_MDSPAN_NAMESPACE::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>;
public:
    using extents_type     = base_type_::extents_type;
    using index_type       = base_type_::index_type;
    using reference        = base_type_::reference;

    using base_type_::accessor;
    using base_type_::data_handle;
    using base_type_::mapping;
    using base_type_::mdspan;
    using base_type_::rank;

    template<class... SizeTypes , typename = std::enable_if_t<(( 
        rank() == sizeof...(SizeTypes) &&
        (detail::valid_indices<index_type, SizeTypes...>()) ))>>
        __attribute__((always_inline)) constexpr reference operator()(SizeTypes... indices) const {
        return accessor().access(data_handle(), mapping()(static_cast<index_type>(std::move(indices))...));
    }

    template<class SizeType , typename = std::enable_if_t<((
        std::is_convertible_v<const SizeType&, index_type> &&
        std::is_nothrow_constructible_v<index_type, const SizeType&> ))>>
        __attribute__((always_inline)) constexpr reference operator()(const std::array<SizeType, rank()>& indices) const {
        return operator[](indices);
    }
};

// Deduction guide
template<class ElementType, class... SizeTypes, typename = std::enable_if_t<(
            ((std::is_convertible_v<SizeTypes, size_t>) && ...) &&
            (sizeof...(SizeTypes) > 0))>>
HIC_HOST_DEVICE explicit mdspan(ElementType*, SizeTypes...)
  -> mdspan<ElementType, dextents<size_t, sizeof...(SizeTypes)>>;
