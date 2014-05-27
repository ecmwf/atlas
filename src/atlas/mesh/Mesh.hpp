/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef Mesh_hpp
#define Mesh_hpp

#include <map>
#include <vector>
#include <string>

#include "eckit/memory/SharedPtr.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

class FunctionSpace;

//------------------------------------------------------------------------------------------------------

class Mesh {

public: // types

    typedef eckit::SharedPtr<Mesh> Ptr;

public: // methods

    virtual ~Mesh();

    /// checks if function space exists
    bool has_function_space(const std::string& name);

    /// Takes ownership, and will be deleted automatically
    FunctionSpace& add_function_space( FunctionSpace* function_space );

    /// accessor by name
    FunctionSpace& function_space(const std::string& name);

    /// accessor by index
    FunctionSpace& function_space(int idx);

    int nb_function_spaces() { return function_spaces_.size(); }

private: // members

    std::map< std::string, size_t > index_; ///< index of function spaces

    std::vector< FunctionSpace* > function_spaces_; ///< function spaces

};

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" 
{
  Mesh* atlas__Mesh__new ();
  void atlas__Mesh__delete (Mesh* This);
  void atlas__Mesh__add_function_space (Mesh* This, FunctionSpace* function_space); 
  FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name); 
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // Mesh_hpp
