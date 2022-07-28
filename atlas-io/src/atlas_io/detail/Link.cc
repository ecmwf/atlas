/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Link.h"

#include "eckit/filesystem/PathName.h"

#include "atlas_io/RecordItem.h"

namespace atlas {
namespace io {

//bool Link::relative() const {
//    std::string path = RecordItem::URI{uri}.path;
//    if( path.size() == 0 ) {
//        return true;
//    }
//    if( path[0] == '/' ) {
//        return false;
//    }
//    if( path[0] == '~' ) {
//        return false;
//    }
//    return true;
//}

//Link operator/( const eckit::PathName& dir, const Link& link ) {
//    auto relative_path = eckit::PathName{ RecordItem::URI{link}.path };
//    auto absolute_uri = RecordItem::URI(link);
//    absolute_uri.path = dir / relative_path;
//    return Link{absolute_uri.str()};
//};


}  // namespace io
}  // namespace atlas
