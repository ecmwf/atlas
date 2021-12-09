/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas_trace.h"

#include <string>
#include <vector>

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/trace/CodeLocation.h"

namespace atlas {

//static std::vector<std::string> Labels( int size, const char* labels[] ) {
//    std::vector<std::string> _labels;
//    _labels.reserve( size );
//    for ( int i = 0; i < size; ++i ) {
//        _labels.emplace_back( labels[i] );
//    }
//    return _labels;
//}

extern "C" {

Trace* new_atlas_Trace(const char* file, int line, const char* title) {
    return new Trace(CodeLocation(file, line, nullptr, true), std::string(title));
}

Trace* new_atlas_Trace_labels_1(const char* file, int line, const char* title, const char* label1) {
    std::vector<std::string> labels{label1};
    return new Trace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}
Trace* new_atlas_Trace_labels_2(const char* file, int line, const char* title, const char* label1, const char* label2) {
    std::vector<std::string> labels{label1, label2};
    return new Trace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}
Trace* new_atlas_Trace_labels_3(const char* file, int line, const char* title, const char* label1, const char* label2,
                                const char* label3) {
    std::vector<std::string> labels{label1, label2, label3};
    return new Trace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}
Trace* new_atlas_Trace_labels_4(const char* file, int line, const char* title, const char* label1, const char* label2,
                                const char* label3, const char* label4) {
    std::vector<std::string> labels{label1, label2, label3, label4};
    return new Trace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}
Trace* new_atlas_Trace_labels_5(const char* file, int line, const char* title, const char* label1, const char* label2,
                                const char* label3, const char* label4, const char* label5) {
    std::vector<std::string> labels{label1, label2, label3, label4, label5};
    return new Trace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}

void delete_atlas_Trace(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot delete uninitialised atlas_Trace");
    delete This;
    This = nullptr;
}

void atlas_Trace__start(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot start uninitialised atlas_Trace");
    This->start();
}

void atlas_Trace__stop(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot stop uninitialised atlas_Trace");
    This->stop();
}

void atlas_Trace__pause(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot pause uninitialised atlas_Trace");
    This->pause();
}

void atlas_Trace__resume(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot resume uninitialised atlas_Trace");
    This->resume();
}

int atlas_Trace__running(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot check 'running' status of uninitialised atlas_Trace");
    return This->running();
}

double atlas_Trace__elapsed(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot check elapsed time of uninitialised atlas_Trace");
    return This->elapsed();
}

}  // extern C

}  // namespace atlas
