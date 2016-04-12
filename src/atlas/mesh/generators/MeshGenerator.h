/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/util/Config.h"

namespace atlas {
namespace mesh {
    class Mesh;
} }

namespace atlas {
namespace grid {
    class Grid;
    class GridDistribution;
} }

namespace atlas {
namespace mesh {
namespace generators {

// -----------------------------------------------------------------------------

class MeshGenerator : public eckit::Owned {

public:

  typedef eckit::SharedPtr<MeshGenerator> Ptr;
  typedef atlas::util::Config Parameters;
  static MeshGenerator* create(const std::string &, const eckit::Parametrisation & = Parameters());

public:

    MeshGenerator();

    virtual ~MeshGenerator();

    virtual void generate( const grid::Grid&, const grid::GridDistribution&, Mesh& ) const =0;
    virtual void generate( const grid::Grid&, Mesh& ) const =0;

    Mesh* generate( const grid::Grid&, const grid::GridDistribution& ) const;
    Mesh* generate( const grid::Grid& ) const;

    Mesh* operator()( const grid::Grid&, const grid::GridDistribution& ) const;
    Mesh* operator()( const grid::Grid& ) const;

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

// -----------------------------------------------------------------------------

#define Parametrisation eckit::Parametrisation
#define grid_Grid grid::Grid
#define grid_GridDistribution grid::GridDistribution

extern "C" {
void atlas__MeshGenerator__delete(MeshGenerator* This);
MeshGenerator* atlas__MeshGenerator__create(const char* name, const Parametrisation* params);
Mesh* atlas__MeshGenerator__generate__grid_griddist(const MeshGenerator* This, const grid_Grid* grid, const grid_GridDistribution* distribution);
Mesh* atlas__MeshGenerator__generate__grid(const MeshGenerator* This, const grid_Grid* grid);
}

#undef grid_Grid
#undef grid_GridDistribution
#undef Parametrisation

} // namespace generators
} // namespace mesh
} // namespace atlas

#endif
