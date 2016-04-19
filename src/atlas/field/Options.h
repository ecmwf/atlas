/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_field_Options_h
#define atlas_field_Options_h

#include "atlas/internals/atlas_config.h"
#include "atlas/util/Config.h"

// ----------------------------------------------------------------------------
// Forward declarations

// ----------------------------------------------------------------------------

namespace atlas {
namespace field {

// ----------------------------------------------------------------------------

class global : public util::Config
{
public:
  global(size_t _owner=0)
  {
    set("global",true);
    set("owner",_owner);
  }
};

class levels : public util::Config
{
public:
  levels(size_t _levels)
  {
    set("levels",_levels);
  }
};

class name : public util::Config
{
public:
  name(const std::string& _name)
  {
    set("name",_name);
  }
};

template< typename T >
class datatypeT : public util::Config
{
public:
  datatypeT()
  {
    set("datatype",array::DataType::kind<T>());
  }
};


class datatype : public util::Config
{
public:
  datatype(array::DataType::kind_t kind)
  {
    set("datatype",kind);
  }
  datatype(const std::string& str)
  {
    set("datatype",array::DataType::str_to_kind(str));
  }

};

// ----------------------------------------------------------------------------

} // namespace field
} // namespace atlas

#endif // atlas_field_Options_h
