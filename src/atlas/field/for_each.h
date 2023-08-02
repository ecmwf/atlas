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

#include "atlas/field/detail/for_each.h"

#include "atlas/array/Array.h"

#include "Field.h"
#include "FieldSet.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace field {

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_value_masked( const eckit::Parametrisation& , const Mask& , std::tuple<Field...>&& , const Function& )

template <typename Mask, typename... Field, typename Function>
void for_each_value_masked(const eckit::Parametrisation& config, const Mask& mask, std::tuple<Field...>&& fields, const Function& function) {
    auto field_1 = std::get<0>(fields);
    if constexpr (std::is_same_v<std::decay_t<Mask>,atlas::Field>) {
        auto h_dim = field_1.horizontal_dimension();

        ATLAS_ASSERT( mask.datatype() == array::make_datatype<int>() );
        ATLAS_ASSERT( mask.rank() <= h_dim.size() );

        if (h_dim.size() == 1) {
            ATLAS_ASSERT(h_dim[0] == 0);
            auto mask_view = array::make_view<const int,1>(mask);
            auto mask_wrap = [mask_view](idx_t i, auto&&... args) { return mask_view(i); };
            return for_each_value_masked(config, mask_wrap, std::move(fields), function);
        }
        else if (h_dim.size() == 2) {
            if (mask.rank() == 1) {
                auto mask_view_1d      = array::make_view<const int,1>(mask);
                auto mask_view_shape2d = array::make_shape(field_1.shape(h_dim[0]), field_1.shape(h_dim[1]));
                auto mask_view         = array::View<const int,2>( mask_view_1d.data(), mask_view_shape2d );
                if( h_dim[0] == 0 && h_dim[1] == 2) {
                    auto mask_wrap = [mask_view](idx_t i, idx_t /*dummy*/, idx_t j, auto&&... args) { return mask_view(i,j); };
                    return for_each_value_masked(config, mask_wrap, std::move(fields), function);
                }
                else {
                    ATLAS_ASSERT(h_dim[0] == 0 && h_dim[1] == 1);
                    auto mask_wrap = [mask_view](idx_t i, idx_t j, auto&&... args) { return mask_view(i,j); };
                    return for_each_value_masked(config, mask_wrap, std::move(fields), function);
                }
            }
            else {
                auto mask_view = array::make_view<const int,2>(mask);
                if( h_dim[0] == 0 && h_dim[1] == 2) {
                    auto mask_wrap = [mask_view](idx_t i, idx_t /*dummy*/, idx_t j, auto&&... args) { return mask_view(i,j); };
                    return for_each_value_masked(config, mask_wrap, std::move(fields), function);
                }
                else {
                    ATLAS_ASSERT(h_dim[0] == 0 && h_dim[1] == 1);
                    auto mask_wrap = [mask_view](idx_t i, idx_t j, auto&&... args) { return mask_view(i,j); };
                    return for_each_value_masked(config, mask_wrap, std::move(fields), function);
                }
            }
        }
        else {
            ATLAS_THROW_EXCEPTION("More than 2 horizontal indices is not yet supported");
        }
    }
    else {
        constexpr auto num_fields = std::tuple_size_v<std::tuple<Field...>>;
        static_assert(num_fields == detail::function_traits<Function>::arity,"!");
        using value_type = std::decay_t<detail::first_argument<Function>>;
        switch (field_1.rank()) {
            case 1: return detail::for_each_value_masked_rank<value_type,1>(config,mask,std::move(fields),function);
            case 2: return detail::for_each_value_masked_rank<value_type,2>(config,mask,std::move(fields),function);
            case 3: return detail::for_each_value_masked_rank<value_type,3>(config,mask,std::move(fields),function);
            case 4: return detail::for_each_value_masked_rank<value_type,4>(config,mask,std::move(fields),function);
            case 5: return detail::for_each_value_masked_rank<value_type,5>(config,mask,std::move(fields),function);
            default: ATLAS_THROW_EXCEPTION("Only fields with rank <= 5 are currently supported. Given rank: " << std::get<0>(fields).rank());
        }
    }
    ATLAS_THROW_EXCEPTION("Invalid function");
}

template <typename Mask, typename Function>
void for_each_value_masked(const eckit::Parametrisation& config, const Mask& mask, Field field, const Function& function) {
    return for_each_value_masked(config, mask, std::make_tuple(field), function);
}

template <typename Mask, typename Function>
void for_each_value_masked(const eckit::Parametrisation& config, const Mask& mask, Field field_1, Field field_2, const Function& function) {
    return for_each_value_masked(config, mask, std::make_tuple(field_1, field_2), function);
}

template <typename Mask, typename Function>
void for_each_value_masked(const eckit::Parametrisation& config, const Mask& mask, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_value_masked(config, mask, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename Mask, typename Function>
void for_each_value_masked(const eckit::Parametrisation& config, const Mask& mask, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_value_masked(config, mask, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_value_masked( const ExecutionPolicy&& , const Mask& , std::tuple<Field...>&& , const Function& )

template <typename ExecutionPolicy, typename Mask, typename... Field, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value_masked(ExecutionPolicy, const Mask& mask, std::tuple<Field...>&& fields, const Function& function) {
    return for_each_value_masked(option::execution_policy<ExecutionPolicy>(), mask, std::move(fields), function);
}

template <typename ExecutionPolicy, typename Mask, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value_masked(ExecutionPolicy execution_policy, const Mask& mask, Field field, const Function& function) {
    return for_each_value_masked(execution_policy, mask, std::make_tuple(field), function);
}

template <typename ExecutionPolicy, typename Mask, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value_masked(ExecutionPolicy execution_policy, const Mask& mask, Field field_1, Field field_2, const Function& function) {
    return for_each_value_masked(execution_policy, mask, std::make_tuple(field_1, field_2), function);
}

template <typename ExecutionPolicy, typename Mask, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value_masked(ExecutionPolicy execution_policy, const Mask& mask, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_value_masked(execution_policy, mask, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename ExecutionPolicy, typename Mask, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value_masked(ExecutionPolicy execution_policy, const Mask& mask, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_value_masked(execution_policy, mask, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_value_masked( const Mask& , std::tuple<Field...>&& , const Function& )

template <typename Mask, typename... Field, typename Function>
void for_each_value_masked(const Mask& mask, std::tuple<Field...>&& fields, const Function& function) {
    return for_each_value_masked(util::NoConfig(), mask, std::move(fields), function);
}

template <typename Mask, typename Function>
void for_each_value_masked(const Mask& mask, Field field, const Function& function) {
    return for_each_value_masked(mask, std::make_tuple(field), function);
}

template <typename Mask, typename Function>
void for_each_value_masked(const Mask& mask, Field field_1, Field field_2, const Function& function) {
    return for_each_value_masked(mask, std::make_tuple(field_1, field_2), function);
}

template <typename Mask, typename Function>
void for_each_value_masked(const Mask& mask, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_value_masked(mask, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename Mask, typename Function>
void for_each_value_masked(const Mask& mask, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_value_masked(mask, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_value( const eckit::Parametrisation& , std::tuple<Field...>&& , const Function& )

template <typename... Field, typename Function>
void for_each_value(const eckit::Parametrisation& config, std::tuple<Field...>&& fields, const Function& function) {
    return for_each_value_masked(config, array::helpers::detail::no_mask, std::move(fields), function);
}

template <typename Function>
void for_each_value(const eckit::Parametrisation& config, Field field, const Function& function) {
    return for_each_value_masked(config, std::make_tuple(field), function);
}

template <typename Function>
void for_each_value(const eckit::Parametrisation& config, Field field_1, Field field_2, const Function& function) {
    return for_each_value_masked(config, std::make_tuple(field_1, field_2), function);
}

template <typename Function>
void for_each_value(const eckit::Parametrisation& config, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_value_masked(config, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename Function>
void for_each_value(const eckit::Parametrisation& config, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_value_masked(config, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_value( const ExecutionPolicy&& , std::tuple<Field...>&& , const Function& )

template <typename ExecutionPolicy, typename... Field, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value(ExecutionPolicy, std::tuple<Field...>&& fields, const Function& function) {
    return for_each_value_masked(option::execution_policy<ExecutionPolicy>(), array::helpers::detail::no_mask, std::move(fields), function);
}

template <typename ExecutionPolicy, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value(ExecutionPolicy execution_policy, Field field, const Function& function) {
    return for_each_value(execution_policy, std::make_tuple(field), function);
}

template <typename ExecutionPolicy, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value(ExecutionPolicy execution_policy, Field field_1, Field field_2, const Function& function) {
    return for_each_value(execution_policy, std::make_tuple(field_1, field_2), function);
}

template <typename ExecutionPolicy, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value(ExecutionPolicy execution_policy, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_value(execution_policy, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename ExecutionPolicy, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_value(ExecutionPolicy execution_policy, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_value(execution_policy, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_value( std::tuple<Field...>&& , const Function& )

template <typename... Field, typename Function>
void for_each_value(std::tuple<Field...>&& fields, const Function& function) {
    return for_each_value_masked(util::NoConfig(), array::helpers::detail::no_mask, std::move(fields), function);
}

template <typename Function>
void for_each_value(Field field, const Function& function) {
    return for_each_value(std::make_tuple(field), function);
}

template <typename Function>
void for_each_value(Field field_1, Field field_2, const Function& function) {
    return for_each_value(std::make_tuple(field_1, field_2), function);
}

template <typename Function>
void for_each_value(Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_value(std::make_tuple(field_1, field_2, field_3), function);
}

template <typename Function>
void for_each_value(Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_value(std::make_tuple(field_1, field_2, field_3, field_4), function);
}



//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_column_masked( const eckit::Parametrisation& , const Mask& , std::tuple<Field...>&& , const Function& )

template <typename Mask, typename... Field, typename Function>
void for_each_column_masked(const eckit::Parametrisation& config, const Mask& mask, std::tuple<Field...>&& fields, const Function& function) {
    if constexpr (std::is_same_v<std::decay_t<Mask>,atlas::Field>) {
        auto field_1 = std::get<0>(fields);
        auto h_dim = field_1.horizontal_dimension();

        ATLAS_ASSERT( mask.datatype() == array::make_datatype<int>() );
        ATLAS_ASSERT( mask.rank() <= h_dim.size() );

        if (h_dim.size() == 1) {
            return for_each_column_masked(config, array::make_view<const int,1>(mask), std::move(fields), function);
        }
        else if (h_dim.size() == 2) {
            if (mask.rank() == 1) {
                auto mask_view_1d      = array::make_view<const int,1>(mask);
                auto mask_view_shape2d = array::make_shape(field_1.shape(h_dim[0]), field_1.shape(h_dim[1]));
                auto mask_view         = array::View<const int,2>( mask_view_1d.data(), mask_view_shape2d );
                return for_each_column_masked(config, mask_view, std::move(fields), function);
            }
            else {
                return for_each_column_masked(config, array::make_view<const int,2>(mask), std::move(fields), function);
            }
        }
        else {
            ATLAS_THROW_EXCEPTION("More than 2 horizontal indices is not yet supported");
        }
    }
    else {
        switch (std::get<0>(fields).rank()) {
            case 2: return detail::for_each_column_masked_rank<2>(config, mask, std::move(fields), function);
            case 3: return detail::for_each_column_masked_rank<3>(config, mask, std::move(fields), function);
            case 4: return detail::for_each_column_masked_rank<4>(config, mask, std::move(fields), function);
        }
    }
    ATLAS_THROW_EXCEPTION("Invalid function");
}

template <typename Mask, typename Function>
void for_each_column_masked(const eckit::Parametrisation& config, const Mask& mask, Field field, const Function& function) {
    return for_each_column_masked(config, mask, std::make_tuple(field), function);
}

template <typename Mask, typename Function>
void for_each_column_masked(const eckit::Parametrisation& config, const Mask& mask, Field field_1, Field field_2, const Function& function) {
    return for_each_column_masked(config, mask, std::make_tuple(field_1, field_2), function);
}

template <typename Mask, typename Function>
void for_each_column_masked(const eckit::Parametrisation& config, const Mask& mask, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_column_masked(config, mask, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename Mask, typename Function>
void for_each_column_masked(const eckit::Parametrisation& config, const Mask& mask, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_column_masked(config, mask, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_column_masked( const ExecutionPolicy&& , const Mask& , std::tuple<Field...>&& , const Function& )

template <typename ExecutionPolicy, typename Mask, typename... Field, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column_masked(ExecutionPolicy, const Mask& mask, std::tuple<Field...>&& fields, const Function& function) {
    return for_each_column_masked(option::execution_policy<ExecutionPolicy>(), mask, std::move(fields), function);
}

template <typename ExecutionPolicy, typename Mask, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column_masked(ExecutionPolicy execution_policy, const Mask& mask, Field field, const Function& function) {
    return for_each_column_masked(execution_policy, mask, std::make_tuple(field), function);
}

template <typename ExecutionPolicy, typename Mask, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column_masked(ExecutionPolicy execution_policy, const Mask& mask, Field field_1, Field field_2, const Function& function) {
    return for_each_column_masked(execution_policy, mask, std::make_tuple(field_1, field_2), function);
}

template <typename ExecutionPolicy, typename Mask, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column_masked(ExecutionPolicy execution_policy, const Mask& mask, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_column_masked(execution_policy, mask, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename ExecutionPolicy, typename Mask, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column_masked(ExecutionPolicy execution_policy, const Mask& mask, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_column_masked(execution_policy, mask, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_column_masked( const Mask& , std::tuple<Field...>&& , const Function& )

template <typename Mask, typename... Field, typename Function>
void for_each_column_masked(const Mask& mask, std::tuple<Field...>&& fields, const Function& function) {
    return for_each_column_masked(util::NoConfig(), mask, std::move(fields), function);
}

template <typename Mask, typename Function>
void for_each_column_masked(const Mask& mask, Field field, const Function& function) {
    return for_each_column_masked(mask, std::make_tuple(field), function);
}

template <typename Mask, typename Function>
void for_each_column_masked(const Mask& mask, Field field_1, Field field_2, const Function& function) {
    return for_each_column_masked(mask, std::make_tuple(field_1, field_2), function);
}

template <typename Mask, typename Function>
void for_each_column_masked(const Mask& mask, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_column_masked(mask, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename Mask, typename Function>
void for_each_column_masked(const Mask& mask, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_column_masked(mask, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_column( const eckit::Parametrisation& , std::tuple<Field...>&& , const Function& )

template <typename... Field, typename Function>
void for_each_column(const eckit::Parametrisation& config, std::tuple<Field...>&& fields, const Function& function) {
    return for_each_column_masked(config, array::helpers::detail::no_mask, std::move(fields), function);
}

template <typename Function>
void for_each_column(const eckit::Parametrisation& config, Field field, const Function& function) {
    return for_each_column_masked(config, std::make_tuple(field), function);
}

template <typename Function>
void for_each_column(const eckit::Parametrisation& config, Field field_1, Field field_2, const Function& function) {
    return for_each_column_masked(config, std::make_tuple(field_1, field_2), function);
}

template <typename Function>
void for_each_column(const eckit::Parametrisation& config, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_column_masked(config, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename Function>
void for_each_column(const eckit::Parametrisation& config, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_column_masked(config, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_column( const ExecutionPolicy&& , std::tuple<Field...>&& , const Function& )

template <typename ExecutionPolicy, typename... Field, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column(ExecutionPolicy, std::tuple<Field...>&& fields, const Function& function) {
    return for_each_column_masked(option::execution_policy<ExecutionPolicy>(), array::helpers::detail::no_mask, std::move(fields), function);
}

template <typename ExecutionPolicy, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column(ExecutionPolicy execution_policy, Field field, const Function& function) {
    return for_each_column(execution_policy, std::make_tuple(field), function);
}

template <typename ExecutionPolicy, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column(ExecutionPolicy execution_policy, Field field_1, Field field_2, const Function& function) {
    return for_each_column(execution_policy, std::make_tuple(field_1, field_2), function);
}

template <typename ExecutionPolicy, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column(ExecutionPolicy execution_policy, Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_column(execution_policy, std::make_tuple(field_1, field_2, field_3), function);
}

template <typename ExecutionPolicy, typename Function, typename = std::enable_if_t<execution::is_execution_policy<ExecutionPolicy>()>>
void for_each_column(ExecutionPolicy execution_policy, Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_column(execution_policy, std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------
// void for_each_column( std::tuple<Field...>&& , const Function& )

template <typename... Field, typename Function>
void for_each_column(std::tuple<Field...>&& fields, const Function& function) {
    return for_each_column_masked(util::NoConfig(), array::helpers::detail::no_mask, std::move(fields), function);
}

template <typename Function>
void for_each_column(Field field, const Function& function) {
    return for_each_column(std::make_tuple(field), function);
}

template <typename Function>
void for_each_column(Field field_1, Field field_2, const Function& function) {
    return for_each_column(std::make_tuple(field_1, field_2), function);
}

template <typename Function>
void for_each_column(Field field_1, Field field_2, Field field_3, const Function& function) {
    return for_each_column(std::make_tuple(field_1, field_2, field_3), function);
}

template <typename Function>
void for_each_column(Field field_1, Field field_2, Field field_3, Field field_4, const Function& function) {
    return for_each_column(std::make_tuple(field_1, field_2, field_3, field_4), function);
}

//------------------------------------------------------------------------------------------------------------------------------------------------

} // namespace field
} // namespace atlas
