#pragma once

// This file contains preprocessor black magic. It contains macros that
// can return the number of arguments passed

//-----------------------------------------------------------------------------------------------------------
// Public

/// Returns the number of passed arguments
#define __ATLAS__NARG(...)

/// Splice a and b together
#define __ATLAS__SPLICE(a,b)

//-----------------------------------------------------------------------------------------------------------
// Details

// Undefine these, to be redefined further down.
#undef __ATLAS__NARG
#undef __ATLAS__SPLICE

#define __ATLAS__REVERSE 5, 4, 3, 2, 1, 0
#define __ATLAS__ARGN(_1, _2, _3, _4, _5, N, ...) N
#define __ATLAS__NARG__(dummy, ...)  __ATLAS__ARGN(__VA_ARGS__)
#define __ATLAS__NARG_(...)          __ATLAS__NARG__(dummy, ##__VA_ARGS__, __ATLAS__REVERSE)
#define __ATLAS__SPLICE(a,b)         __ATLAS__SPLICE_1(a,b)
#define __ATLAS__SPLICE_1(a,b)       __ATLAS__SPLICE_2(a,b)
#define __ATLAS__SPLICE_2(a,b)       a##b


#define __ATLAS__ARG16(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, ...) _15
#define __ATLAS__HAS_COMMA(...) __ATLAS__ARG16(__VA_ARGS__, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
#define __ATLAS__TRIGGER_PARENTHESIS(...) ,
#define __ATLAS__ISEMPTY(...)                                                                      \
__ATLAS__ISEMPTY_(                                                                                 \
  /* test if there is just one argument, eventually an empty one */                                \
  __ATLAS__HAS_COMMA(__VA_ARGS__),                                                                 \
  /* test if _TRIGGER_PARENTHESIS_ together with the argument adds a comma */                      \
  __ATLAS__HAS_COMMA(__ATLAS__TRIGGER_PARENTHESIS __VA_ARGS__),                                    \
  /* test if the argument together with a parenthesis adds a comma */                              \
  __ATLAS__HAS_COMMA(__VA_ARGS__ (/*empty*/)),                                                     \
  /* test if placing it between __ATLAS__TRIGGER_PARENTHESIS and the parenthesis adds a comma */   \
  __ATLAS__HAS_COMMA(__ATLAS__TRIGGER_PARENTHESIS __VA_ARGS__ (/*empty*/))                         \
)

#define __ATLAS__PASTE5(_0, _1, _2, _3, _4) _0 ## _1 ## _2 ## _3 ## _4
#define __ATLAS__ISEMPTY_(_0, _1, _2, _3) __ATLAS__HAS_COMMA(__ATLAS__PASTE5(__ATLAS__IS_EMPTY_CASE_, _0, _1, _2, _3))
#define __ATLAS__IS_EMPTY_CASE_0001 ,

#define __ATLAS__CALL_NARG_1(...) 0
#define __ATLAS__CALL_NARG_0 __ATLAS__NARG_
#define __ATLAS__NARG(...) __ATLAS__SPLICE( __ATLAS__CALL_NARG_, __ATLAS__ISEMPTY( __VA_ARGS__ ) ) (__VA_ARGS__)
