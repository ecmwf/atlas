/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas_io/detail/Time.h"
#include "atlas_io/detail/Version.h"

namespace atlas {
namespace io {

struct RecordInfo {
    Version version_;
    Time created_;

    const Version& version() const { return version_; }
    const Time& created() const { return created_; }
};

}  // namespace io
}  // namespace atlas
