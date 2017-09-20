/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/util/Config.h"
#include "atlas/util/Earth.h"
#include "atlas/array/DataType.h"

// ----------------------------------------------------------------------------

namespace atlas {
namespace option {

// ----------------------------------------------------------------------------

class type : public util::Config
{
public:
  type( const std::string& _type )
  {
    set("type",_type);
  }
};

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

class variables : public util::Config
{
public:
  variables(size_t _variables)
  {
    set("variables",_variables);
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
  datatype(array::DataType dtype)
  {
    set("datatype",dtype.kind());
  }
};

class halo : public util::Config
{
public:
  halo(size_t size)
  {
    set("halo",size);
  }
};

class radius : public util::Config
{
public:
  radius(double _radius)
  {
    set("radius",_radius);
  }
  radius(const std::string& key = "Earth")
  {
    if( key == "Earth" ) {
      set("radius",util::Earth::radiusInMeters());
    } else {
      NOTIMP;
    }
  }

};


// ----------------------------------------------------------------------------

} // namespace option
} // namespace atlas
