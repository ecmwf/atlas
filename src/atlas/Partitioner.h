/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_Partitioner_h
#define atlas_Partitioner_h

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas {

class Grid;
class GridDistribution;

class Partitioner : public eckit::Owned {
public:

  typedef eckit::SharedPtr<Partitioner> Ptr;

public:

  Partitioner(const Grid& grid);
  Partitioner(const Grid& grid, const size_t nb_partitions);
  virtual ~Partitioner();

  virtual void partition( int part[] ) const = 0;

  virtual GridDistribution* distribution() const;

public:

  size_t nb_partitions() const;
  const Grid& grid() const { return grid_; }

private:

  size_t nb_partitions_;
  const Grid& grid_;
};


// ------------------------------------------------------------------

class PartitionerFactory {
  public:
    /*!
     * \brief build Partitioner with factory key, constructor arguments
     * \return Partitioner
     */
    static Partitioner* build(const std::string&, const Grid& grid);
    static Partitioner* build(const std::string&, const Grid& grid, const size_t nb_partitions);

    /*!
     * \brief list all registered field creators
     */
    static void list(std::ostream &);
    static bool has(const std::string& name);

  private:
    std::string name_;
    virtual Partitioner* make(const Grid& grid) = 0 ;
    virtual Partitioner* make(const Grid& grid, const size_t nb_partitions) = 0 ;
  protected:

    PartitionerFactory(const std::string&);
    virtual ~PartitionerFactory();
};

// ------------------------------------------------------------------

template<class T>
class PartitionerBuilder : public PartitionerFactory {
  virtual Partitioner* make(const Grid& grid) {
      return new T(grid);
  }
  virtual Partitioner* make(const Grid& grid, const size_t nb_partitions) {
      return new T(grid, nb_partitions);
  }
  public:
    PartitionerBuilder(const std::string& name) : PartitionerFactory(name) {}
};

// ------------------------------------------------------------------

} // namespace atlas

#endif // atlas_Partitioner_h
