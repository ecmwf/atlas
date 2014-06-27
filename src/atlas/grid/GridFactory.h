/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grid_GridFactory_H
#define atlas_grid_GridFactory_H

#include <map>
#include <string>

#include "atlas/grid/Grid.h"
#include "atlas/grid/GridSpec.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
/// Factory class, that allows Grid's to created from GridSpecs
/// The Grid class derivatives, must register a creator
/// This relies on automatic type registration. This is not well defined for
/// all compilers, and hence derived Grid class, may need to manually
/// register with the factory.
///
/// Requires that the Grid derivatives provide a default constructor.
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8

/// Abstract base class that can create a Grid
/// The constructor will register itself with the GridFactory
class GridCreator {
public:
   GridCreator(const std::string& grid_type);
   virtual ~GridCreator();

   virtual Grid::Ptr create() const = 0;
};


/// Concrete templatized class for creating Grid derivatives
/// These should be registered with GridFactory
template <class T>
class CreatorImpl : public GridCreator {
public:
    CreatorImpl(const std::string& grid_type) : GridCreator(grid_type) {}

    /// Grid derivatives must provide a default constructor
    virtual Grid::Ptr create() const { return Grid::Ptr(new T); }
};


/// GridFactory: responsible for creating Grid derivatives from a GridSpec
class GridFactory {
public:
   /// Will find creator given a GridSpec. Will ask creator to create the GRID.
   /// Once created, the Grid is initialised from the GridSpec and returned.
   static Grid::Ptr create(const GridSpec&);

   /// Register the Creator for a given grid_type.
   /// The derived classes must remember to register with this factory
   static void registerit(const std::string& grid_type, GridCreator* creator);

private:
   static std::map<std::string, GridCreator*>& get_table();
};

/// Place in the header of the Derived Grid class
#define REGISTER(classname) \
private: \
   static const CreatorImpl<classname> creator;

/// Place in source file of the derived Grid class
#define REGISTERIMPL(classname,grid_type) \
   const CreatorImpl<classname> classname::creator(grid_type);

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
