/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_numerics_Nabla_h
#define atlas_numerics_Nabla_h

#include "eckit/memory/Owned.h"

namespace eckit { class Parametrisation; }
namespace atlas { namespace next { class FunctionSpace; } }
namespace atlas { class Field; }

namespace atlas {
namespace numerics {

class Nabla : public eckit::Owned {
public:

  static Nabla* create(const next::FunctionSpace &);
  static Nabla* create(const next::FunctionSpace &, const eckit::Parametrisation &);

public:
  Nabla(const next::FunctionSpace &, const eckit::Parametrisation &);
  virtual ~Nabla();

  virtual void gradient(const Field &scalar, Field &grad) const = 0;
  virtual void divergence(const Field &vector, Field &div) const = 0;
  virtual void curl(const Field &vector, Field &curl) const = 0;
  virtual void laplacian(const Field &scalar, Field &laplacian) const = 0;

private:
  const next::FunctionSpace& fs_;
  const eckit::Parametrisation& config_;

};


// ------------------------------------------------------------------

class NablaFactory {
public:
    /*!
     * \brief build Nabla with factory key, constructor arguments
     * \return Nabla
     */
    static Nabla* build(const next::FunctionSpace &, const eckit::Parametrisation &);

    /*!
     * \brief list all registered field creators
     */
    static void list(std::ostream &);
    static bool has(const std::string& name);

private:
    virtual Nabla* make(const next::FunctionSpace &, const eckit::Parametrisation &) = 0 ;

protected:
    NablaFactory(const std::string&);
    virtual ~NablaFactory();

private:
  std::string name_;
};

// ------------------------------------------------------------------

template<class T>
class NablaBuilder : public NablaFactory {

public:
    NablaBuilder(const std::string& name) : NablaFactory(name) {}

private:
    virtual Nabla* make(const next::FunctionSpace & fs, const eckit::Parametrisation & p) {
        return new T(fs,p);
    }
};

// ------------------------------------------------------------------
#define NextFunctionSpace next::FunctionSpace
#define Parametrisation eckit::Parametrisation
extern "C" {

void atlas__Nabla__delete (Nabla* This);
Nabla* atlas__Nabla__create (const NextFunctionSpace* functionspace, const Parametrisation* params);
void atlas__Nabla__gradient (const Nabla* This, const Field* scalar, Field* grad);
void atlas__Nabla__divergence (const Nabla* This, const Field* vector, Field* div);
void atlas__Nabla__curl (const Nabla* This, const Field* vector, Field* curl);
void atlas__Nabla__laplacian (const Nabla* This, const Field* scalar, Field* laplacian);
}
#undef NextFunctionSpace
#undef Parametrisation

} // namespace numerics
} // namespace atlas

#endif // atlas_numerics_Nabla_h
