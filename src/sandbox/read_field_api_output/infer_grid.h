/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstddef>
#include <string>
#include <map>

 namespace atlas {
 inline std::string infer_grid_from_nproma_and_nblks(int nproma, int nblks) {
    std::map<size_t,std::string> size_to_grid = {
        {2048,  "F16"},
        {6114,  "N32"},
        {35718, "N80"},
    };
    std::string grid;
    for (const auto & [k,v] : size_to_grid) {
        if (nproma * nblks >= k && nproma * (nblks-1) < k) {
            grid = v;
            break;
        }
    }
    return grid;
 }
 }
