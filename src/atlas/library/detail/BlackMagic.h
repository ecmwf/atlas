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

// This file contains preprocessor black magic. It contains macros that
// can return the number of arguments passed

//-----------------------------------------------------------------------------------------------------------
// Public

/// Returns the number of passed arguments
#define __ATLAS_NARG(...)

/// Splice a and b together
#define __ATLAS_SPLICE(a, b)

#define __ATLAS_STRINGIFY(a) a

#define __ATLAS_TYPE(Type, ...)
#define __ATLAS_TYPE_SCOPE(Type, ...)

//-----------------------------------------------------------------------------------------------------------
// Details

// Undefine these, to be redefined further down.
#undef __ATLAS_NARG
#undef __ATLAS_SPLICE
#undef __ATLAS_TYPE
#undef __ATLAS_TYPE_SCOPE

#define __ATLAS_REVERSE 5, 4, 3, 2, 1, 0
#define __ATLAS_ARGN(_1, _2, _3, _4, _5, N, ...) N
#define __ATLAS_NARG__(dummy, ...) __ATLAS_ARGN(__VA_ARGS__)
#define __ATLAS_NARG_(...) __ATLAS_NARG__(dummy, ##__VA_ARGS__, __ATLAS_REVERSE)
#define __ATLAS_SPLICE(a, b) __ATLAS_SPLICE_1(a, b)
#define __ATLAS_SPLICE_1(a, b) __ATLAS_SPLICE_2(a, b)
#define __ATLAS_SPLICE_2(a, b) a##b

#define __ATLAS_ARG16(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, ...) _15
#define __ATLAS_HAS_COMMA(...) __ATLAS_ARG16(__VA_ARGS__, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
#define __ATLAS_TRIGGER_PARENTHESIS(...) ,
#define __ATLAS_ISEMPTY(...)                                                                                       \
    __ATLAS_ISEMPTY_(                                /* test if there is just one argument, eventually an empty  \
                      one */ \
                     __ATLAS_HAS_COMMA(__VA_ARGS__), /* test if                  \
                                                      _TRIGGER_PARENTHESIS_    \
                                                      together with the        \
                                                      argument adds a comma */                                 \
                     __ATLAS_HAS_COMMA(                                                                            \
                         __ATLAS_TRIGGER_PARENTHESIS __VA_ARGS__), /* test if the argument together with \
                                            a parenthesis adds a comma */         \
                     __ATLAS_HAS_COMMA(__VA_ARGS__(/*empty*/)),    /* test if placing it between              \
                                       __ATLAS_TRIGGER_PARENTHESIS and the     \
                                       parenthesis adds a comma */    \
                     __ATLAS_HAS_COMMA(__ATLAS_TRIGGER_PARENTHESIS __VA_ARGS__(/*empty*/)))

#define __ATLAS_PASTE5(_0, _1, _2, _3, _4) _0##_1##_2##_3##_4
#define __ATLAS_ISEMPTY_(_0, _1, _2, _3) __ATLAS_HAS_COMMA(__ATLAS_PASTE5(__ATLAS_IS_EMPTY_CASE_, _0, _1, _2, _3))
#define __ATLAS_IS_EMPTY_CASE_0001 ,

#define __ATLAS_NARG(...) __ATLAS_SPLICE(__ATLAS_CALL_NARG_, __ATLAS_ISEMPTY(__VA_ARGS__))(__VA_ARGS__)
#define __ATLAS_CALL_NARG_1(...) 0
#define __ATLAS_CALL_NARG_0 __ATLAS_NARG_

#define __ATLAS_COMMA_ARGS(...) __ATLAS_SPLICE(__ATLAS_COMMA_ARGS_, __ATLAS_ISEMPTY(__VA_ARGS__))(__VA_ARGS__)
#define __ATLAS_COMMA_ARGS_1(...)
#define __ATLAS_COMMA_ARGS_0(...) , __VA_ARGS__

#define __ATLAS_ARGS_OR_DUMMY(...)                                       \
    __ATLAS_SPLICE(__ATLAS_ARGS_OR_DUMMY_, __ATLAS_ISEMPTY(__VA_ARGS__)) \
    (__VA_ARGS__)
#define __ATLAS_ARGS_OR_DUMMY_0(...) __VA_ARGS__
#define __ATLAS_ARGS_OR_DUMMY_1(...) 0

#define __ATLAS_TYPE(Type, ...)                                 \
    __ATLAS_SPLICE(__ATLAS_TYPE_, __ATLAS_ISEMPTY(__VA_ARGS__)) \
    (Type, __ATLAS_ARGS_OR_DUMMY(__VA_ARGS__))
#define __ATLAS_TYPE_1(Type, dummy) Type __ATLAS_SPLICE(__variable_, __LINE__)
#define __ATLAS_TYPE_0(Type, ...) Type __ATLAS_SPLICE(__variable_, __LINE__)(__VA_ARGS__)

#define __ATLAS_TYPE_SCOPE(Type, ...)                                 \
    __ATLAS_SPLICE(__ATLAS_TYPE_SCOPE_, __ATLAS_ISEMPTY(__VA_ARGS__)) \
    (Type, __ATLAS_ARGS_OR_DUMMY(__VA_ARGS__))
#define __ATLAS_TYPE_SCOPE_1(Type, ...)                                                              \
    for (bool __ATLAS_SPLICE(__done_, __LINE__) = false; __ATLAS_SPLICE(__done_, __LINE__) != true;) \
        for (Type __ATLAS_SPLICE(__variable_, __LINE__); __ATLAS_SPLICE(__done_, __LINE__) != true;  \
             __ATLAS_SPLICE(__done_, __LINE__) = true)
#define __ATLAS_TYPE_SCOPE_0(Type, ...)                                                                          \
    for (bool __ATLAS_SPLICE(__done_, __LINE__) = false; __ATLAS_SPLICE(__done_, __LINE__) != true;)             \
        for (Type __ATLAS_SPLICE(__variable_, __LINE__)(__VA_ARGS__); __ATLAS_SPLICE(__done_, __LINE__) != true; \
             __ATLAS_SPLICE(__done_, __LINE__) = true)
