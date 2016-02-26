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
namespace atlas { namespace numerics { class Method; } }
namespace atlas { namespace field { class Field; } }

namespace atlas {
namespace numerics {

class Nabla : public eckit::Owned {
public:

  static Nabla* create(const Method &);
  static Nabla* create(const Method &, const eckit::Parametrisation &);

public:
  Nabla(const Method &, const eckit::Parametrisation &);
  virtual ~Nabla();

  virtual void gradient(const field::Field &scalar, field::Field &grad) const = 0;
  virtual void divergence(const field::Field &vector, field::Field &div) const = 0;
  virtual void curl(const field::Field &vector, field::Field &curl) const = 0;
  virtual void laplacian(const field::Field &scalar, field::Field &laplacian) const = 0;

};


// ------------------------------------------------------------------

class NablaFactory {
public:
    /*!
     * \brief build Nabla with factory key, constructor arguments
     * \return Nabla
     */
    static Nabla* build(const Method &, const eckit::Parametrisation &);

    /*!
     * \brief list all registered field creators
     */
    static void list(std::ostream &);
    static bool has(const std::string& name);

private:
    virtual Nabla* make(const Method &, const eckit::Parametrisation &) = 0 ;

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
    virtual Nabla* make(const Method &method, const eckit::Parametrisation &p) {
        return new T(method,p);
    }
};

// ------------------------------------------------------------------
#define Parametrisation eckit::Parametrisation
#define field_Field field::Field
extern "C" {

void atlas__Nabla__delete (Nabla* This);
Nabla* atlas__Nabla__create (const Method* method, const Parametrisation* params);
void atlas__Nabla__gradient (const Nabla* This, const field_Field* scalar, field_Field* grad);
void atlas__Nabla__divergence (const Nabla* This, const field_Field* vector, field_Field* div);
void atlas__Nabla__curl (const Nabla* This, const field_Field* vector, field_Field* curl);
void atlas__Nabla__laplacian (const Nabla* This, const field_Field* scalar, field_Field* laplacian);
}
#undef field_Field
#undef Parametrisation

} // namespace numerics
} // namespace atlas

#endif // atlas_numerics_Nabla_h
