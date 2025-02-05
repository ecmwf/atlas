/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <type_traits>
#include <utility>
#include <sstream>

#include "atlas/field/Field.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace field {
namespace detail {

template <typename Function> // primary template
struct function_traits : 
    public function_traits<decltype(&std::remove_reference<Function>::type::operator())> { };

template <typename ClassType, typename ReturnType, typename... Arguments>
struct function_traits<ReturnType(ClassType::*)(Arguments...) const> : 
    function_traits<ReturnType(*)(Arguments...)> { };

template <typename ClassType, typename ReturnType, typename... Arguments>
struct function_traits<ReturnType(ClassType::*)(Arguments...)> : 
    function_traits<ReturnType(*)(Arguments...)> { };

template <typename ReturnType, typename... Arguments>
struct function_traits<ReturnType(*)(Arguments...)> {
    using result_type = ReturnType;

    template <std::size_t Index>
    using arg_t = typename std::tuple_element<
        Index,
        std::tuple<Arguments...>
    >::type;

    static constexpr std::size_t arity = sizeof...(Arguments);
};

template <typename Function>
using first_argument = typename function_traits<Function>::template arg_t<0>;

template <typename Function, typename Value>
constexpr bool valid_value_function() {
    return std::is_same_v<std::decay_t<first_argument<Function>>,std::decay_t<Value>>;
}

template <typename Function>
using function_argument_data_type = typename std::decay_t<typename std::decay_t<first_argument<Function>>::value_type>;

template <typename Function>
constexpr idx_t function_argument_rank() {
    return std::decay_t<first_argument<Function>>::rank();
}

template <typename Function, typename Value, int Rank>
constexpr bool valid_column_function() {
    using value_type = std::decay_t<Value>;
    return function_argument_rank<Function>() == Rank && std::is_same_v<function_argument_data_type<Function>,value_type>;
}

template <typename Config,  typename Mask, typename Function, typename... Views, std::size_t... Dims>
void for_each_value_masked_view(std::index_sequence<Dims...>, const Config& config, const Mask& mask, const Function& function, std::tuple<Views...>&& views ) {
    if constexpr (std::is_invocable_r_v<int, Mask, decltype (Dims)...>) {
        array::helpers::ArrayForEach<Dims...>::apply(config,std::move(views),mask,function);
    }
    else {
        ATLAS_THROW_EXCEPTION("Invalid mask function passed");
    }
}

template <int FieldIdx, typename Value, int Rank, typename... Field>
auto make_view_tuple(std::tuple<Field...>&& fields) {
    constexpr auto num_fields = std::tuple_size_v<std::tuple<Field...>>;
    if constexpr (FieldIdx < num_fields) {
        return std::tuple_cat(std::make_tuple(
            array::make_view<Value, Rank>(std::get<FieldIdx>(fields))),
            make_view_tuple<FieldIdx + 1, Value, Rank>(std::move(fields)));
    } else {
        return std::make_tuple();
    }
    ATLAS_UNREACHABLE();
}

template <typename Value, int Rank, typename Config,  typename Mask, typename... Field, typename Function>
void for_each_value_masked_rank(const Config& config, const Mask& mask, std::tuple<Field...>&& fields, const Function& function ) {
    constexpr auto dims = std::make_index_sequence<Rank>();
    if constexpr (valid_value_function<Function,Value>()) {
        return for_each_value_masked_view(dims, config, mask, function, make_view_tuple<0, Value, Rank>(std::move(fields)));
    }
}

template <typename Config,  typename Mask, typename Function, typename... Views >
void for_each_column_masked_view(const Config& config, const Mask& mask, const std::vector<idx_t>& h_dim, const Function& function, std::tuple<Views...>&& views ) {
    using View = std::decay_t<std::tuple_element_t<0, std::tuple<Views...>>>;
    using value_type = typename View::value_type;
    constexpr int rank = View::rank();

    if (h_dim.size()>=rank) {
        ATLAS_THROW_EXCEPTION("Cannot extract column for Rank="<<rank<<" field and " << h_dim.size() << " horizontal indices");
    }
    for (auto h : h_dim) {
        if (h < 0 || h >= rank) {
            ATLAS_THROW_EXCEPTION("Invalid horizontal_dimension " << h_dim[0] << " must less rank " << rank);
        }
    }

    auto has_duplicates = [](const auto &v) {
        for (const auto& i : v) {
            for (const auto & j : v) {
                if (&i == &j) break;
                if (i == j) return true;
            }
        }
        return false;
    };
    if (has_duplicates(h_dim)) {
        ATLAS_THROW_EXCEPTION("horizontal_dimension contains duplicates");
    }

    if constexpr (rank==1) {
        ATLAS_THROW_EXCEPTION("Cannot use for_each_column with Rank=1 fields");
    }
    if constexpr (rank==2) {
        if constexpr (valid_column_function<Function,value_type,rank-1>() &&
                      std::is_invocable_r_v<int, Mask, idx_t>) {
            switch (h_dim[0]) {
                case 0: array::helpers::ArrayForEach<0>::apply(config,std::move(views),mask,function); return;
                case 1: array::helpers::ArrayForEach<1>::apply(config,std::move(views),mask,function); return;
                default: break;
            }
            ATLAS_THROW_EXCEPTION("Not implemented for horizontal_dimension = " << h_dim);
        }
        ATLAS_THROW_EXCEPTION("Invalid function passed");
    }
    if constexpr (rank==3) {
        if (h_dim.size() == 1) {
            if constexpr (valid_column_function<Function,value_type,rank-1>() &&
                          std::is_invocable_r_v<int, Mask, idx_t>) {
                switch (h_dim[0]) {
                    case 0: array::helpers::ArrayForEach<0>::apply(config,std::move(views),mask,function); return;
                    case 1: array::helpers::ArrayForEach<1>::apply(config,std::move(views),mask,function); return;
                    case 2: array::helpers::ArrayForEach<2>::apply(config,std::move(views),mask,function); return;
                    default: break;
                }
                ATLAS_THROW_EXCEPTION("Not implemented for horizontal_dimension = " << h_dim);
            }
            ATLAS_THROW_EXCEPTION("Invalid function passed");
        }
        else if (h_dim.size() == 2) {
            if constexpr (valid_column_function<Function,value_type,rank-2>() &&
                          std::is_invocable_r_v<int, Mask, idx_t, idx_t>) {
                switch (h_dim[0]) {
                    case 0: {
                        switch (h_dim[1]) {
                            case 1: array::helpers::ArrayForEach<0,1>::apply(config,std::move(views),mask,function); return;
                            case 2: array::helpers::ArrayForEach<0,2>::apply(config,std::move(views),mask,function); return;
                            default: break;
                        }
                        break;
                    }
                    case 1: {
                        switch (h_dim[1]) {
                            case 2: array::helpers::ArrayForEach<1,2>::apply(config,std::move(views),mask,function); return;
                            default: break;
                        }
                        break;
                    }
                }
                ATLAS_THROW_EXCEPTION("Not implemented for horizontal_dimension = " << h_dim);
            }
            ATLAS_THROW_EXCEPTION("Invalid function passed");
        }
        ATLAS_THROW_EXCEPTION("Not implemented for h_dim = " << h_dim);
    }
    if constexpr (rank==4) {
        if (h_dim.size() == 1) {
            if constexpr (valid_column_function<Function,value_type,rank-1>() &&
                          std::is_invocable_r_v<int, Mask, idx_t>) {
                switch (h_dim[0]) {
                    case 0: array::helpers::ArrayForEach<0>::apply(config,std::move(views),mask,function); return;
                    case 1: array::helpers::ArrayForEach<1>::apply(config,std::move(views),mask,function); return;
                    case 2: array::helpers::ArrayForEach<2>::apply(config,std::move(views),mask,function); return;
                    case 3: array::helpers::ArrayForEach<3>::apply(config,std::move(views),mask,function); return;
                    default: break;
                }
            }
        }
        else if (h_dim.size() == 2) {
            if constexpr (valid_column_function<Function,value_type,rank-2>() &&
                          std::is_invocable_r_v<int, Mask, idx_t, idx_t>) {
                switch (h_dim[0]) {
                    case 0: {
                        switch (h_dim[1]) {
                            case 1: array::helpers::ArrayForEach<0,1>::apply(config,std::move(views),mask,function); return;
                            case 2: array::helpers::ArrayForEach<0,2>::apply(config,std::move(views),mask,function); return;
                            case 3: array::helpers::ArrayForEach<0,3>::apply(config,std::move(views),mask,function); return;
                            default: break;
                        }
                        break;
                    }
                    case 1: {
                        switch (h_dim[1]) {
                            case 2: array::helpers::ArrayForEach<1,2>::apply(config,std::move(views),mask,function); return;
                            case 3: array::helpers::ArrayForEach<1,3>::apply(config,std::move(views),mask,function); return;
                            default: break;
                        }
                        break;
                    }
                    case 2: {
                        switch (h_dim[1]) {
                            case 3: array::helpers::ArrayForEach<2,3>::apply(config,std::move(views),mask,function); return;
                            default: break;
                        }
                        break;
                    }
                    default: break;
                }
                ATLAS_THROW_EXCEPTION("Not implemented for horizontal_dimension = " << h_dim);
            }
            ATLAS_THROW_EXCEPTION("Invalid function passed");
        }
        ATLAS_THROW_EXCEPTION("Not implemented for h_dim = " << h_dim);
    }
    ATLAS_THROW_EXCEPTION("Not implemented for rank="<<rank);
}

template <int Rank, typename Mask, typename... Fields, typename Function>
void for_each_column_masked_rank(const eckit::Parametrisation& config, const Mask& mask, std::tuple<Fields...>&& fields, const Function& function) {

    constexpr auto num_fields = std::tuple_size_v<std::tuple<Fields...>>;
    ATLAS_ASSERT(num_fields == detail::function_traits<Function>::arity);

    using value_type = detail::function_argument_data_type<Function>;
    ATLAS_ASSERT( std::get<0>(fields).datatype() == array::make_datatype<value_type>() );

    auto h_dim = std::get<0>(fields).horizontal_dimension();
    ATLAS_ASSERT(Rank ==  detail::function_argument_rank<Function>() + h_dim.size());

    return for_each_column_masked_view(config, mask, h_dim, function, make_view_tuple<0, value_type, Rank>(std::move(fields)));
}

} // namespace detail
} // namespace field
} // namespace atlas
