/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/method/Method.h"

#include <map>
#include "eckit/exception/Exceptions.h"
#include "eckit/linalg/LinearAlgebra.h"
#include "eckit/linalg/Vector.h"
#include "eckit/log/Timer.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/thread/Once.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace interpolation {

namespace {


typedef std::map<std::string, MethodFactory*> MethodFactoryMap_t;
static MethodFactoryMap_t *m = 0;
static eckit::Mutex *local_mutex = 0;
static pthread_once_t once = PTHREAD_ONCE_INIT;


static void init() {
    local_mutex = new eckit::Mutex();
    m = new MethodFactoryMap_t();
}


}  // (anonymous namespace)


MethodFactory::MethodFactory(const std::string& name):
    name_(name) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    if (m->find(name) != m->end()) {
        throw eckit::SeriousBug("MethodFactory duplicate '" + name + "'");
    }

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


MethodFactory::~MethodFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


Method* MethodFactory::build(const std::string& name, const Method::Config& config) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    MethodFactoryMap_t::const_iterator j = m->find(name);
    if (j == m->end()) {
        eckit::Log::error() << "MethodFactory '" << name << "' not found." << std::endl;
        eckit::Log::error() << "MethodFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j) {
            eckit::Log::error() << '\t' << (*j).first << std::endl;
        }
        throw eckit::SeriousBug("MethodFactory '" + name + "' not found.");
    }

    return (*j).second->make(config);
}


void Method::execute(const FieldSet& fieldsSource, FieldSet& fieldsTarget) const {
    eckit::TraceTimer<Atlas> tim("atlas::interpolation::method::Method::execute()");

    const size_t N = fieldsSource.size();
    ASSERT(N == fieldsTarget.size());

    for (size_t i = 0; i < fieldsSource.size(); ++i) {
        Log::debug<Atlas>() << "Method::execute() on field " << (i+1) << '/' << N << "..." << std::endl;

        const Field& src = fieldsSource[i];
        Field& tgt = fieldsTarget[i];

        eckit::linalg::Vector v_src( const_cast<double*>(src.data<double>()), src.shape(0) );
        eckit::linalg::Vector v_tgt( tgt.data<double>(), tgt.shape(0) );

        eckit::linalg::LinearAlgebra::backend().spmv(matrix_, v_src, v_tgt);
    }
}


void Method::execute(const Field& fieldSource, Field& fieldTarget) const {
    eckit::TraceTimer<Atlas> tim("atlas::interpolation::method::Method::execute()");

    eckit::linalg::Vector
            v_src(const_cast< Field& >(fieldSource).data<double>(), fieldSource.shape(0)),
            v_tgt(fieldTarget.data<double>(), fieldTarget.shape(0));

    eckit::linalg::LinearAlgebra::backend().spmv(matrix_, v_src, v_tgt);
}


void Method::normalise(Triplets& triplets) {
    // sum all calculated weights for normalisation
    double sum = 0.0;

    for (size_t j = 0; j < triplets.size(); ++j) {
        sum += triplets[j].value();
    }

    // now normalise all weights according to the total
    const double invSum = 1.0 / sum;
    for (size_t j = 0; j < triplets.size(); ++j) {
        triplets[j].value() *= invSum;
    }
}

}  // interpolation
}  // atlas

