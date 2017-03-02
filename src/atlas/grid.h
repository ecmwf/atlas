/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck

#pragma once

#include <map>
#include <string>

#include "atlas/grid/Grid.h"

namespace atlas {
namespace grid {

class GridCreator {

public:

  using Registry = std::map< std::string, GridCreator*>;
  static const Registry& nameRegistry();
  static const Registry& typeRegistry();

public:

  GridCreator( const std::string& name );
  GridCreator( const std::string& type, const std::string& name );

  ~GridCreator();

  virtual const Grid::grid_t* create( const Grid::Config& ) const;

  virtual const Grid::grid_t* create( const std::string& ) const =0;

  std::string name() const;

  std::string type() const;

protected:

  bool match(const std::string& string, std::vector<std::string>& matches) const;

private:

  friend std::ostream& operator<<( std::ostream& os, const GridCreator& g );

  virtual void print(std::ostream& os) const =0;

  std::string name_;
  std::string type_;

};

}  // namespace grid
}  // namespace atlas
