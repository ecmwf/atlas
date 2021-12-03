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

// C wrapper interfaces to C++ routines
extern "C" {
void atlas__atlas_init_noargs();
void atlas__atlas_finalize();
const char* atlas__eckit_version();
const char* atlas__eckit_git_sha1();
const char* atlas__eckit_git_sha1_abbrev(int length);
const char* atlas__atlas_version();
const char* atlas__atlas_git_sha1();
const char* atlas__atlas_git_sha1_abbrev(int length);
}
