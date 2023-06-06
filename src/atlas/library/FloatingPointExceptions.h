/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <csignal>
#include <iosfwd>
#include <map>
#include <string>
#include <cfenv>

namespace atlas {
namespace library {

// ------------------------------------------------------------------------------------

/* @brief Enable floating point exceptions
 * 
 * An environment variable is responsible to enable:
 * - ATLAS_FPE=0                  --> Default. do not enable anything in Atlas (other libraries or programs may still do so)
 * - ATLAS_FPE=1                  --> FE_INVALID, FE_DIVBYZERO, FE_OVERFLOW
 * - ATLAS_FPE=VALUE1,VALUE2,...  --> Specified codes only. Valid codes: 
 *                                    FE_INVALID, FE_DIVBYZERO, FE_OVERFLOW, FE_UNDERFLOW, FE_INEXACT, FE_ALL_EXCEPT
 * 
 * @note This function is called automatically within Library::initialize() and should
 *       not be called directly
 */
void enable_floating_point_exceptions();

/* @brief Enable floating point exception
 * 
 * Valid codes: 
 *   FE_INVALID, FE_DIVBYZERO, FE_OVERFLOW, FE_UNDERFLOW, FE_INEXACT, FE_ALL_EXCEPT
 * 
 * @return false when the exception was already enabled, true when change was made
 */
bool enable_floating_point_exception(const std::string&);
bool enable_floating_point_exception(int);

/* @brief Disable floating point exception
 * 
 * Valid codes: 
 *   FE_INVALID, FE_DIVBYZERO, FE_OVERFLOW, FE_UNDERFLOW, FE_INEXACT, FE_ALL_EXCEPT
 * 
 * @return false when the exception was already disabled, true when change was made
 */
bool disable_floating_point_exception(const std::string&);
bool disable_floating_point_exception(int);

// ------------------------------------------------------------------------------------

/* @brief Enable atlas signal handler for all signals
 * 
 * An environment variable is responsible to enable:
 * - ATLAS_SIGNAL_HANDLER=0   --> Default. Atlas will not set any signal handlers
 * - ATLAS_SIGNAL_HANDLER=1   --> Enable atlas_signal_handler instead of the default for signals:
 * 
 * Enabled signals:
 * - SIGABRT
 * - SIGFPE
 * - SIGILL
 * - SIGINT
 * - SIGSEGV
 * - SIGTERM
 * - SIGKILL
 * 
 * @note This function is called automatically within Library::initialize() and should
 *       not be called directly
 */
void enable_atlas_signal_handler();

// ------------------------------------------------------------------------------------

using signal_handler_t = void (*)(int);
using signal_action_t  = void (*)(int, ::siginfo_t*, void*);

[[noreturn]] void atlas_signal_handler(int signum, ::siginfo_t* si, void* unused);

// ------------------------------------------------------------------------------------

}  // namespace library
}  // namespace atlas
