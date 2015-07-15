/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_meshgen_MeshGenerator_h
#define atlas_meshgen_MeshGenerator_h

#include <iosfwd>
#include <string>

#include "eckit/memory/Owned.h"
#include "eckit/config/Parametrisation.h"
#include "atlas/Config.h"

namespace atlas {

class Mesh;
class Grid;
class GridDistribution;

namespace meshgen {

//------------------------------------------------------------------------------------------------------

class MeshGenerator : public eckit::Owned {

public:

    typedef atlas::Config Parameters; // temporary until Parameters class exists

    MeshGenerator();

    virtual ~MeshGenerator();

    virtual void generate( const Grid&, const GridDistribution&, Mesh& ) const =0;
    virtual void generate( const Grid&, Mesh& ) const =0;

    Mesh* operator()( const Grid&, const GridDistribution& ) const;
    Mesh* operator()( const Grid& ) const;

};



class MeshGeneratorFactory {
  public:
    /*!
     * \brief build MeshGenerator with factory key, and default options
     * \return mesh generator
     */
    static MeshGenerator* build(const std::string&);

    /*!
     * \brief build MeshGenerator with factory key inside parametrisation,
     * and options specified in parametrisation as well
     * \return mesh generator
     */
    static MeshGenerator* build(const std::string&, const eckit::Parametrisation&);

    /*!
     * \brief list all registered mesh generators
     */
    static void list(std::ostream &);

  private:
    std::string name_;
    virtual MeshGenerator* make() = 0 ;
    virtual MeshGenerator* make(const eckit::Parametrisation&) = 0 ;

  protected:

    MeshGeneratorFactory(const std::string&);
    virtual ~MeshGeneratorFactory();


};


template<class T>
class MeshGeneratorBuilder : public MeshGeneratorFactory {
  virtual MeshGenerator* make() {
      return new T();
  }
  virtual MeshGenerator* make(const eckit::Parametrisation& param) {
        return new T(param);
    }
  public:
    MeshGeneratorBuilder(const std::string& name) : MeshGeneratorFactory(name) {}
};

//------------------------------------------------------------------------------------------------------

} // namespace meshgen
} // namespace atlas

#endif
