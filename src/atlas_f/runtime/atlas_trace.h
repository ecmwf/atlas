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

namespace atlas {
class Trace;
}

extern "C" {

atlas::Trace* new_atlas_Trace(const char* file, int line, const char* title);
atlas::Trace* new_atlas_Trace_labels_1(const char* file, int line, const char* title, const char* label1);
atlas::Trace* new_atlas_Trace_labels_2(const char* file, int line, const char* title, const char* label1,
                                       const char* label2);
atlas::Trace* new_atlas_Trace_labels_3(const char* file, int line, const char* title, const char* label1,
                                       const char* label2, const char* label3);
atlas::Trace* new_atlas_Trace_labels_4(const char* file, int line, const char* title, const char* label1,
                                       const char* label2, const char* label3, const char* label4);
atlas::Trace* new_atlas_Trace_labels_5(const char* file, int line, const char* title, const char* label1,
                                       const char* label2, const char* label3, const char* label4, const char* label5);
void delete_atlas_Trace(atlas::Trace* This);
void atlas_Trace__start(atlas::Trace* This);
void atlas_Trace__stop(atlas::Trace* This);
void atlas_Trace__pause(atlas::Trace* This);
void atlas_Trace__resume(atlas::Trace* This);
int atlas_Trace__running(atlas::Trace* This);
double atlas_Trace__elapsed(atlas::Trace* This);

}  // extern C
