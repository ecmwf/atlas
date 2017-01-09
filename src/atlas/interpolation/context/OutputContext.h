/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_interpolation_context_OutputContext_h
#define atlas_interpolation_context_OutputContext_h

#include "atlas/interpolation/context/Context.h"


namespace atlas {
namespace interpolation {
namespace context {


struct OutputContext : Context {

    OutputContext(
            const std::string& gridname,
            const std::string& partitioner,
            const std::string& meshGenerator,
            const mesh::Mesh::Ptr prePartitionedMesh,
            const grid::Domain& prePartitionedDomain = grid::Domain::makeGlobal(),
            bool meshGeneratorTriangulate = false,
            double meshGeneratorAngle = 0 );

    virtual void write(const std::string& fileName) /*= 0;*/ { NOTIMP; }

};


struct OutputContextFactory {
    static OutputContext* build(const std::string& key, const std::string& name);
protected:
    std::string name_;
    OutputContextFactory(const std::string&);
    virtual ~OutputContextFactory();
    virtual OutputContext* make(const std::string& name) = 0;
};


template<class T>
struct OutputContextBuilder : public OutputContextFactory {
    OutputContextBuilder(const std::string& name) : OutputContextFactory(name) {}
private:
    virtual OutputContext* make(const std::string& name) { return new T(name); }
};


}  // context
}  // interpolation
}  // atlas


#endif
