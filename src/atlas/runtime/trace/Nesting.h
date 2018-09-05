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

#include "atlas/runtime/trace/CodeLocation.h"
#include "atlas/runtime/trace/CallStack.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

class Nesting {
public:
    Nesting( const CodeLocation&, const std::string& id = "" );
    ~Nesting();
    operator CallStack() const { return stack_; }
    void stop();
    void start();

private:
    CallStack stack_;
    CodeLocation loc_;
    std::string id_;
    bool running_{true};
};

}  // namespace trace
}  // namespace runtime
}  // namespace atlas