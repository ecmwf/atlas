/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_interpolation_context_InputContext_h
#define atlas_interpolation_context_InputContext_h

#include "atlas/interpolation/context/Context.h"


namespace atlas {
namespace interpolation {
namespace context {


struct InputContext : Context {

    InputContext(
            const std::string& gridname,
            const std::string& partitioner,
            const std::string& meshGenerator,
            bool meshGeneratorTriangulate = false,
            double meshGeneratorAngle = 0 );

    virtual void read(const std::string& name);

};


struct InputContextFactory {
    static InputContext* build(const std::string& key, const std::string& name);
protected:
    std::string name_;
    InputContextFactory(const std::string&);
    virtual ~InputContextFactory();
    virtual InputContext* make(const std::string& name) = 0;
};


template<class T>
struct InputContextBuilder : public InputContextFactory {
    InputContextBuilder(const std::string& name) : InputContextFactory(name) {}
private:
    virtual InputContext* make(const std::string& name) { return new T(name); }
};


}  // context
}  // interpolation
}  // atlas


#endif
