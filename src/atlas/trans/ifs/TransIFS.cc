/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/log/JSON.h"

#include "atlas/library/config.h"
#if ATLAS_HAVE_ECTRANS
#include "ectrans/transi.h"
#else
#include "transi/trans.h"
#endif

#include "atlas/array.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/detail/TransFactory.h"
#include "atlas/trans/ifs/TransIFS.h"

#include "atlas/grid/detail/partitioner/CheckerboardPartitioner.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"

using Topology = atlas::mesh::Nodes::Topology;
using atlas::Field;
using atlas::FunctionSpace;
using atlas::array::ArrayView;
using atlas::array::LocalView;
using atlas::array::make_shape;
using atlas::array::make_view;
using atlas::functionspace::NodeColumns;
using atlas::functionspace::Spectral;
using atlas::functionspace::StructuredColumns;
using atlas::mesh::IsGhostNode;

namespace atlas {
namespace trans {

namespace {
static TransBuilderGrid<TransIFS> builder_ifs("ifs", "ifs");  // deprecated
static TransBuilderGrid<TransIFS> builder_ectrans("ectrans", "ectrans");
}  // namespace

class TransParameters {
public:
    TransParameters(const TransIFS& trans, const eckit::Configuration& config): trans_(trans), config_(config) {}
    ~TransParameters() = default;

    bool scalar_derivatives() const { return config_.getBool("scalar_derivatives", false); }

    bool wind_EW_derivatives() const { return config_.getBool("wind_EW_derivatives", false); }

    bool vorticity_divergence_fields() const { return config_.getBool("vorticity_divergence_fields", false); }

    bool split_y() const { return config_.getBool("split_y", true); }

    int fft() const {
        static const std::map<std::string, int> string_to_FFT = {{"FFT992", TRANS_FFT992}, {"FFTW", TRANS_FFTW}};
        return string_to_FFT.at(config_.getString("fft", "FFTW"));
    }

    bool flt() const { return config_.getBool("flt", false); }

    std::string read_legendre() const { return config_.getString("read_legendre", ""); }

    std::string write_legendre() const { return config_.getString("write_legendre", ""); }

    bool global() const { return config_.getBool("global", false); }

    int nproma() const { return config_.getInt("nproma", trans_->ngptot); }

    int ngpblks() const {
        int _ngptot = trans_->ngptot;
        int _nproma = nproma();
        ATLAS_ASSERT(_ngptot % _nproma == 0);  // assert _ngptot is divisable by nproma
        return _ngptot / _nproma;
    }

private:
    const Trans_t* trans_;
    const eckit::Configuration& config_;
};

namespace {
std::string fieldset_functionspace(const FieldSet& fields) {
    std::string functionspace("undefined");
    for (idx_t jfld = 0; jfld < fields.size(); ++jfld) {
        if (functionspace == "undefined") {
            functionspace = fields[jfld].functionspace().type();
        }
        if (fields[jfld].functionspace().type() != functionspace) {
            throw_Exception(": fielset has fields with different functionspaces", Here());
        }
    }
    return functionspace;
}
void assert_spectral_functionspace(const FieldSet& fields) {
    for (idx_t jfld = 0; jfld < fields.size(); ++jfld) {
        ATLAS_ASSERT(functionspace::Spectral(fields[jfld].functionspace()));
    }
}

void trans_check(const int code, const char* msg, const eckit::CodeLocation& location) {
    if (code != TRANS_SUCCESS) {
        std::stringstream errmsg;
        errmsg << "atlas::trans ERROR: " << msg << " failed: \n";
        errmsg << ::trans_error_msg(code);
        throw_Exception(errmsg.str(), location);
    }
}
#define TRANS_CHECK(CALL) trans_check(CALL, #CALL, Here())


static int compute_nfld(const Field& f) {
    int nfld = 1;
    for (int i = 1; i < f.rank(); ++i) {
        nfld *= f.shape(i);
    }
    return nfld;
}

static int compute_nfld(const FieldSet& f) {
    int nfld = 0;
    for (int i = 0; i < f.size(); ++i) {
        nfld += compute_nfld(f[i]);
    }
    return nfld;
}

// --------------------------------------------------------------------------------------------

}  // namespace

void TransIFS::dirtrans(const Field& gpfield, Field& spfield, const eckit::Configuration& config) const {
    ATLAS_ASSERT(functionspace::Spectral(spfield.functionspace()));
    if (functionspace::StructuredColumns(gpfield.functionspace())) {
        __dirtrans(functionspace::StructuredColumns(gpfield.functionspace()), gpfield,
                   functionspace::Spectral(spfield.functionspace()), spfield, config);
    }
    else if (functionspace::NodeColumns(gpfield.functionspace())) {
        __dirtrans(functionspace::NodeColumns(gpfield.functionspace()), gpfield,
                   functionspace::Spectral(spfield.functionspace()), spfield, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::dirtrans(const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& config) const {
    ATLAS_TRACE();
    assert_spectral_functionspace(spfields);
    std::string functionspace(fieldset_functionspace(gpfields));

    if (functionspace == StructuredColumns::type()) {
        __dirtrans(StructuredColumns(gpfields[0].functionspace()), gpfields, Spectral(spfields[0].functionspace()),
                   spfields, config);
    }
    else if (functionspace == NodeColumns::type()) {
        __dirtrans(NodeColumns(gpfields[0].functionspace()), gpfields, Spectral(spfields[0].functionspace()), spfields,
                   config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::dirtrans_wind2vordiv(const Field& gpwind, Field& spvor, Field& spdiv,
                                    const eckit::Configuration& config) const {
    ATLAS_ASSERT(Spectral(spvor.functionspace()));
    ATLAS_ASSERT(Spectral(spdiv.functionspace()));
    if (StructuredColumns(gpwind.functionspace())) {
        __dirtrans_wind2vordiv(StructuredColumns(gpwind.functionspace()), gpwind, Spectral(spvor.functionspace()),
                               spvor, spdiv, config);
    }
    else if (NodeColumns(gpwind.functionspace())) {
        __dirtrans_wind2vordiv(NodeColumns(gpwind.functionspace()), gpwind, Spectral(spvor.functionspace()), spvor,
                               spdiv, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}


void TransIFS::dirtrans_adj(const Field& spfield, Field& gpfield, const eckit::Configuration& config) const {
    ATLAS_ASSERT(functionspace::Spectral(spfield.functionspace()));
    if (functionspace::StructuredColumns(gpfield.functionspace())) {
        __dirtrans_adj(functionspace::Spectral(spfield.functionspace()), spfield,
                       functionspace::StructuredColumns(gpfield.functionspace()), gpfield,
                       config);
    }
    else if (functionspace::NodeColumns(gpfield.functionspace())) {
        __dirtrans_adj(functionspace::Spectral(spfield.functionspace()), spfield,
                       functionspace::NodeColumns(gpfield.functionspace()), gpfield,
                       config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::dirtrans_adj(const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& config) const {
    ATLAS_TRACE();
    assert_spectral_functionspace(spfields);
    std::string functionspace(fieldset_functionspace(gpfields));
    if (functionspace == StructuredColumns::type()) {
        __dirtrans_adj(Spectral(spfields[0].functionspace()), spfields,
                       StructuredColumns(gpfields[0].functionspace()), gpfields,
                       config);
    }
    else if (functionspace == NodeColumns::type()) {
 //       __dirtrans_adj(NodeColumns(gpfields[0].functionspace()), gpfields, Spectral(spfields[0].functionspace()), spfields,
 //                 config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::dirtrans_wind2vordiv_adj(const Field& spvor, const Field& spdiv, Field& gpwind,
                                        const eckit::Configuration& config) const {
    ATLAS_ASSERT(Spectral(spvor.functionspace()));
    ATLAS_ASSERT(Spectral(spdiv.functionspace()));
    if (StructuredColumns(gpwind.functionspace())) {
        __dirtrans_wind2vordiv_adj(Spectral(spvor.functionspace()), spvor, spdiv,
                                   StructuredColumns(gpwind.functionspace()), gpwind,
                                   config);
    }
    else if (NodeColumns(gpwind.functionspace())) {
//        __dirtrans_wind2vordiv_adj(Spectral(spvor.functionspace()), spvor, spdiv,
//                                   NodeColumns(gpwind.functionspace()), gpwind,
//                                   config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::invtrans(const Field& spfield, Field& gpfield, const eckit::Configuration& config) const {
    ATLAS_TRACE();
    ATLAS_ASSERT(Spectral(spfield.functionspace()));
    if (StructuredColumns(gpfield.functionspace())) {
        __invtrans(Spectral(spfield.functionspace()), spfield, StructuredColumns(gpfield.functionspace()), gpfield,
                   config);
    }
    else if (NodeColumns(gpfield.functionspace())) {
        __invtrans(Spectral(spfield.functionspace()), spfield, NodeColumns(gpfield.functionspace()), gpfield, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::invtrans(const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& config) const {
    ATLAS_TRACE();
    assert_spectral_functionspace(spfields);
    std::string functionspace(fieldset_functionspace(gpfields));

    if (functionspace == StructuredColumns::type()) {
        __invtrans(Spectral(spfields[0].functionspace()), spfields, StructuredColumns(gpfields[0].functionspace()),
                   gpfields, config);
    }
    else if (functionspace == NodeColumns::type()) {
        __invtrans(Spectral(spfields[0].functionspace()), spfields, NodeColumns(gpfields[0].functionspace()), gpfields,
                   config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}


void TransIFS::invtrans_adj(const Field& gpfield, Field& spfield, const eckit::Configuration& config) const {
    ATLAS_ASSERT(Spectral(spfield.functionspace()));
    if (StructuredColumns(gpfield.functionspace())) {
        __invtrans_adj(Spectral(spfield.functionspace()), spfield, StructuredColumns(gpfield.functionspace()), gpfield,
                       config);
    }
    else if (NodeColumns(gpfield.functionspace())) {
        __invtrans_adj(Spectral(spfield.functionspace()), spfield, NodeColumns(gpfield.functionspace()), gpfield,
                       config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::invtrans_adj(const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& config) const {
    assert_spectral_functionspace(spfields);
    std::string functionspace(fieldset_functionspace(gpfields));

    if (functionspace == StructuredColumns::type()) {
        __invtrans_adj(Spectral(spfields[0].functionspace()), spfields, StructuredColumns(gpfields[0].functionspace()),
                       gpfields, config);
    }
    else if (functionspace == NodeColumns::type()) {
        __invtrans_adj(Spectral(spfields[0].functionspace()), spfields, NodeColumns(gpfields[0].functionspace()),
                       gpfields, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::invtrans_grad(const Field& spfield, Field& gradfield, const eckit::Configuration& config) const {
    ATLAS_ASSERT(Spectral(spfield.functionspace()));
    ATLAS_ASSERT(NodeColumns(gradfield.functionspace()));
    if (StructuredColumns(gradfield.functionspace())) {
        __invtrans_grad(Spectral(spfield.functionspace()), spfield, StructuredColumns(gradfield.functionspace()),
                        gradfield, config);
    }
    else if (NodeColumns(gradfield.functionspace())) {
        __invtrans_grad(Spectral(spfield.functionspace()), spfield, NodeColumns(gradfield.functionspace()), gradfield,
                        config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::invtrans_grad(const FieldSet& spfields, FieldSet& gradfields, const eckit::Configuration& config) const {
    assert_spectral_functionspace(spfields);
    std::string functionspace(fieldset_functionspace(gradfields));

    if (functionspace == StructuredColumns::type()) {
        __invtrans_grad(Spectral(spfields[0].functionspace()), spfields,
                        StructuredColumns(gradfields[0].functionspace()), gradfields, config);
    }
    else if (functionspace == NodeColumns::type()) {
        __invtrans_grad(Spectral(spfields[0].functionspace()), spfields, NodeColumns(gradfields[0].functionspace()),
                        gradfields, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::invtrans_grad_adj(const Field& gradfield, Field& spfield, const eckit::Configuration& config) const {
    ATLAS_ASSERT(Spectral(spfield.functionspace()));
    ATLAS_ASSERT(NodeColumns(gradfield.functionspace()));
    if (StructuredColumns(gradfield.functionspace())) {
        __invtrans_grad_adj(Spectral(spfield.functionspace()), spfield, StructuredColumns(gradfield.functionspace()),
                            gradfield, config);
    }
    else if (NodeColumns(gradfield.functionspace())) {
        __invtrans_grad_adj(Spectral(spfield.functionspace()), spfield, NodeColumns(gradfield.functionspace()),
                            gradfield, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransIFS::invtrans_grad_adj(const FieldSet& gradfields, FieldSet& spfields,
                                 const eckit::Configuration& config) const {
    assert_spectral_functionspace(spfields);
    std::string functionspace(fieldset_functionspace(gradfields));

    if (functionspace == StructuredColumns::type()) {
        __invtrans_grad_adj(Spectral(spfields[0].functionspace()), spfields,
                            StructuredColumns(gradfields[0].functionspace()), gradfields, config);
    }
    else if (functionspace == NodeColumns::type()) {
        __invtrans_grad_adj(Spectral(spfields[0].functionspace()), spfields, NodeColumns(gradfields[0].functionspace()),
                            gradfields, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}



void TransIFS::invtrans_vordiv2wind(const Field& spvor, const Field& spdiv, Field& gpwind,
                                    const eckit::Configuration& config) const {
    ATLAS_ASSERT(Spectral(spvor.functionspace()));
    ATLAS_ASSERT(Spectral(spdiv.functionspace()));
    if (StructuredColumns(gpwind.functionspace())) {
        __invtrans_vordiv2wind(Spectral(spvor.functionspace()), spvor, spdiv, StructuredColumns(gpwind.functionspace()),
                               gpwind, config);
    }
    else if (NodeColumns(gpwind.functionspace())) {
        __invtrans_vordiv2wind(Spectral(spvor.functionspace()), spvor, spdiv, NodeColumns(gpwind.functionspace()),
                               gpwind, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}


void TransIFS::invtrans_vordiv2wind_adj(const Field& gpwind, Field& spvor, Field& spdiv,
                                        const eckit::Configuration& config) const {
    ATLAS_ASSERT(Spectral(spvor.functionspace()));
    ATLAS_ASSERT(Spectral(spdiv.functionspace()));
    if (StructuredColumns(gpwind.functionspace())) {
        __invtrans_vordiv2wind_adj(Spectral(spvor.functionspace()), spvor, spdiv,
                                   StructuredColumns(gpwind.functionspace()), gpwind, config);
    }
    else if (NodeColumns(gpwind.functionspace())) {
        __invtrans_vordiv2wind_adj(Spectral(spvor.functionspace()), spvor, spdiv, NodeColumns(gpwind.functionspace()),
                                   gpwind, config);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::invtrans(const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                        const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                        const eckit::Configuration& config) const {
    ATLAS_TRACE("TransIFS::invtrans");
    TransParameters params(*this, config);
    struct ::InvTrans_t args = new_invtrans(trans_.get());
    args.nscalar             = nb_scalar_fields;
    args.rspscalar           = scalar_spectra;
    args.nvordiv             = nb_vordiv_fields;
    args.rspvor              = vorticity_spectra;
    args.rspdiv              = divergence_spectra;
    args.rgp                 = gp_fields;
    args.lglobal             = params.global();
    args.lscalarders         = params.scalar_derivatives();
    args.luvder_EW           = params.wind_EW_derivatives();
    args.lvordivgp           = params.vorticity_divergence_fields();
    args.nproma              = params.nproma();
    args.ngpblks             = params.ngpblks();
    TRANS_CHECK(::trans_invtrans(&args));
}

void TransIFS::invtrans_adj(const int nb_scalar_fields, const double gp_fields[], const int nb_vordiv_fields,
                            double vorticity_spectra[], double divergence_spectra[], double scalar_spectra[],
                            const eckit::Configuration& config) const {
#if ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ)
    ATLAS_TRACE("TransIFS::invtrans_adj");
    TransParameters params(*this, config);
    struct ::InvTransAdj_t args = new_invtrans_adj(trans_.get());
    args.nscalar                = nb_scalar_fields;
    args.rspscalar              = scalar_spectra;
    args.nvordiv                = nb_vordiv_fields;
    args.rspvor                 = vorticity_spectra;
    args.rspdiv                 = divergence_spectra;
    args.rgp                    = gp_fields;
    args.lglobal                = params.global();
    args.lscalarders            = params.scalar_derivatives();
    args.luvder_EW              = params.wind_EW_derivatives();
    args.lvordivgp              = params.vorticity_divergence_fields();
    args.nproma                 = params.nproma();
    args.ngpblks                = params.ngpblks();
    TRANS_CHECK(::trans_invtrans_adj(&args));
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}


///////////////////////////////////////////////////////////////////////////////

void TransIFS::invtrans(const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                        const eckit::Configuration& config) const {
    return invtrans(nb_scalar_fields, scalar_spectra, 0, nullptr, nullptr, gp_fields, config);
}

///////////////////////////////////////////////////////////////////////////////

void TransIFS::invtrans(const int nb_vordiv_fields, const double vorticity_spectra[], const double divergence_spectra[],
                        double gp_fields[], const eckit::Configuration& config) const {
    return invtrans(0, nullptr, nb_vordiv_fields, vorticity_spectra, divergence_spectra, gp_fields, config);
}

///////////////////////////////////////////////////////////////////////////////


void TransIFS::invtrans_adj(const int nb_scalar_fields, const double gp_fields[], double scalar_spectra[],
                            const eckit::Configuration& config) const {
    return invtrans_adj(nb_scalar_fields, gp_fields, 0, nullptr, nullptr, scalar_spectra, config);
}

///////////////////////////////////////////////////////////////////////////////


void TransIFS::invtrans_adj(const int nb_vordiv_fields, const double gp_fields[], double divergence_spectra[],
                            double vorticity_spectra[], const eckit::Configuration& config) const {
    return invtrans_adj(0, gp_fields, nb_vordiv_fields, vorticity_spectra, divergence_spectra, nullptr, config);
}


void TransIFS::dirtrans(const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                        const eckit::Configuration& config) const {
    ATLAS_TRACE();
    TransParameters params(*this, config);
    struct ::DirTrans_t args = new_dirtrans(trans_.get());
    args.nscalar             = nb_fields;
    args.rgp                 = scalar_fields;
    args.rspscalar           = scalar_spectra;
    args.lglobal             = params.global();
    args.nproma              = params.nproma();
    args.ngpblks             = params.ngpblks();
    TRANS_CHECK(::trans_dirtrans(&args));
}

///////////////////////////////////////////////////////////////////////////////

void TransIFS::dirtrans(const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                        double divergence_spectra[], const eckit::Configuration& config) const {
    ATLAS_TRACE();
    TransParameters params(*this, config);
    struct ::DirTrans_t args = new_dirtrans(trans_.get());
    args.nvordiv             = nb_fields;
    args.rspvor              = vorticity_spectra;
    args.rspdiv              = divergence_spectra;
    args.rgp                 = wind_fields;
    args.lglobal             = params.global();
    args.nproma              = params.nproma();
    args.ngpblks             = params.ngpblks();
    TRANS_CHECK(::trans_dirtrans(&args));
}

}  // namespace trans

// anonymous namespace
namespace {

struct PackNodeColumns {
    LocalView<double, 2>& rgpview_;
    IsGhostNode is_ghost;
    size_t f;

    PackNodeColumns(LocalView<double, 2>& rgpview, const NodeColumns& fs):
        rgpview_(rgpview), is_ghost(fs.nodes()), f(0) {}

    void operator()(const Field& field, idx_t components = 0) {
        switch (field.rank()) {
            case 1:
                pack_1(field, components);
                break;
            case 2:
                pack_2(field, components);
                break;
            case 3:
                pack_3(field, components);
                break;
            default:
                ATLAS_DEBUG_VAR(field.rank());
                ATLAS_NOTIMPLEMENTED;
                //break;
        }
    }

    void pack_1(const Field& field, idx_t) {
        auto gpfield = make_view<double, 1>(field);
        idx_t n      = 0;
        for (idx_t jnode = 0; jnode < gpfield.shape(0); ++jnode) {
            if (!is_ghost(jnode)) {
                rgpview_(f, n) = gpfield(jnode);
                ++n;
            }
        }
        ++f;
    }
    void pack_2(const Field& field, idx_t) {
        auto gpfield      = make_view<double, 2>(field);
        const idx_t nvars = gpfield.shape(1);
        for (idx_t jvar = 0; jvar < nvars; ++jvar) {
            idx_t n = 0;
            for (idx_t jnode = 0; jnode < gpfield.shape(0); ++jnode) {
                if (!is_ghost(jnode)) {
                    rgpview_(f, n) = gpfield(jnode, jvar);
                    ++n;
                }
            }
            ++f;
        }
    }
    void pack_3(const Field& field, idx_t components) {
        auto gpfield = make_view<double, 3>(field);
        if (not components) {
            components = gpfield.shape(2);
        }
        for (idx_t jcomp = 0; jcomp < components; ++jcomp) {
            for (idx_t jlev = 0; jlev < gpfield.shape(1); ++jlev) {
                idx_t n = 0;
                for (idx_t jnode = 0; jnode < gpfield.shape(0); ++jnode) {
                    if (!is_ghost(jnode)) {
                        rgpview_(f, n) = gpfield(jnode, jlev, jcomp);
                        ++n;
                    }
                }
                ++f;
            }
        }
    }
};

struct PackStructuredColumns {
    LocalView<double, 2>& rgpview_;
    size_t f;

    PackStructuredColumns(LocalView<double, 2>& rgpview): rgpview_(rgpview), f(0) {}

    void operator()(const StructuredColumns& sc, const Field& field, idx_t components = 0) {
        switch (field.rank()) {
            case 1:
                pack_1(sc, field, components);
                break;
            case 2:
                pack_2(sc, field, components);
                break;
            case 3:
                pack_3(sc, field, components);
                break;
            default:
                ATLAS_DEBUG_VAR(field.rank());
                ATLAS_NOTIMPLEMENTED;
                //break;
        }
    }

    void pack_1(const StructuredColumns& sc, const Field& field, idx_t) {
        auto gpfield = make_view<double, 1>(field);

        for (idx_t jnode = 0; jnode < sc.sizeOwned(); ++jnode) {
            rgpview_(f, jnode) = gpfield(jnode);
        }
        ++f;
    }
    void pack_2(const StructuredColumns& sc, const Field& field, idx_t) {
        auto gpfield      = make_view<double, 2>(field);
        const idx_t nvars = gpfield.shape(1);
        for (idx_t jvar = 0; jvar < nvars; ++jvar) {
            for (idx_t jnode = 0; jnode < sc.sizeOwned(); ++jnode) {
                rgpview_(f, jnode) = gpfield(jnode, jvar);
            }
            ++f;
        }
    }
    void pack_3(const StructuredColumns& sc, const Field& field, idx_t components) {
        auto gpfield = make_view<double, 3>(field);
        if (not components) {
            components = gpfield.shape(2);
        }
        for (idx_t jcomp = 0; jcomp < components; ++jcomp) {
            const idx_t nvars = gpfield.shape(1);
            for (idx_t jvar = 0; jvar < nvars; ++jvar) {
                for (idx_t jnode = 0; jnode < sc.sizeOwned(); ++jnode) {
                    rgpview_(f, jnode) = gpfield(jnode, jvar, jcomp);
                }
                ++f;
            }
        }
    }
};

struct PackSpectral {
    LocalView<double, 2>& rspecview_;
    size_t f;
    PackSpectral(LocalView<double, 2>& rspecview): rspecview_(rspecview), f(0) {}

    void operator()(const Field& field) {
        switch (field.rank()) {
            case 1:
                pack_1(field);
                break;
            case 2:
                pack_2(field);
                break;
            default:
                ATLAS_DEBUG_VAR(field.rank());
                ATLAS_NOTIMPLEMENTED;
                //break;
        }
    }

    void pack_1(const Field& field) {
        auto spfield = make_view<double, 1>(field);

        for (idx_t jwave = 0; jwave < spfield.shape(0); ++jwave) {
            rspecview_(jwave, f) = spfield(jwave);
        }
        ++f;
    }
    void pack_2(const Field& field) {
        auto spfield = make_view<double, 2>(field);

        const idx_t nvars = spfield.shape(1);

        for (idx_t jvar = 0; jvar < nvars; ++jvar) {
            for (idx_t jwave = 0; jwave < spfield.shape(0); ++jwave) {
                rspecview_(jwave, f) = spfield(jwave, jvar);
            }
            ++f;
        }
    }
};

struct UnpackNodeColumns {
    const LocalView<double, 2>& rgpview_;
    IsGhostNode is_ghost;
    size_t f;

    UnpackNodeColumns(const LocalView<double, 2>& rgpview, const NodeColumns& fs):
        rgpview_(rgpview), is_ghost(fs.nodes()), f(0) {}

    void operator()(Field& field, int components = 0) {
        switch (field.rank()) {
            case 1:
                unpack_1(field, components);
                break;
            case 2:
                unpack_2(field, components);
                break;
            case 3:
                unpack_3(field, components);
                break;
            default:
                ATLAS_DEBUG_VAR(field.rank());
                ATLAS_NOTIMPLEMENTED;
                //break;
        }
    }

    void unpack_1(Field& field, idx_t) {
        auto gpfield = make_view<double, 1>(field);
        idx_t n(0);
        for (idx_t jnode = 0; jnode < gpfield.shape(0); ++jnode) {
            if (!is_ghost(jnode)) {
                gpfield(jnode) = rgpview_(f, n);
                ++n;
            }
        }
        ++f;
    }
    void unpack_2(Field& field, idx_t) {
        auto gpfield      = make_view<double, 2>(field);
        const idx_t nvars = gpfield.shape(1);
        for (idx_t jvar = 0; jvar < nvars; ++jvar) {
            idx_t n = 0;
            for (idx_t jnode = 0; jnode < gpfield.shape(0); ++jnode) {
                if (!is_ghost(jnode)) {
                    gpfield(jnode, jvar) = rgpview_(f, n);
                    ++n;
                }
            }
            ++f;
        }
    }
    void unpack_3(Field& field, idx_t components) {
        auto gpfield = make_view<double, 3>(field);
        if (not components) {
            components = gpfield.shape(2);
        }
        for (idx_t jcomp = 0; jcomp < components; ++jcomp) {
            for (idx_t jlev = 0; jlev < gpfield.shape(1); ++jlev) {
                idx_t n = 0;
                for (idx_t jnode = 0; jnode < gpfield.shape(0); ++jnode) {
                    if (!is_ghost(jnode)) {
                        gpfield(jnode, jlev, jcomp) = rgpview_(f, n);
                        ++n;
                    }
                }
                ++f;
            }
        }
    }
};

struct UnpackStructuredColumns {
    const LocalView<double, 2>& rgpview_;
    size_t f;

    UnpackStructuredColumns(const LocalView<double, 2>& rgpview): rgpview_(rgpview), f(0) {}

    void operator()(const StructuredColumns& sc, Field& field, int components = 0) {
        switch (field.rank()) {
            case 1:
                unpack_1(sc, field, components);
                break;
            case 2:
                unpack_2(sc, field, components);
                break;
            case 3:
                unpack_3(sc, field, components);
                break;
            default:
                ATLAS_DEBUG_VAR(field.rank());
                ATLAS_NOTIMPLEMENTED;
                //break;
        }
    }

    void unpack_1(const StructuredColumns& sc, Field& field, idx_t) {
        auto gpfield = make_view<double, 1>(field);
        for (idx_t jnode = 0; jnode < sc.sizeOwned(); ++jnode) {
            gpfield(jnode) = rgpview_(f, jnode);
        }
        ++f;
    }
    void unpack_2(const StructuredColumns& sc, Field& field, idx_t) {
        auto gpfield      = make_view<double, 2>(field);
        const idx_t nvars = gpfield.shape(1);
        for (idx_t jvar = 0; jvar < nvars; ++jvar) {
            for (idx_t jnode = 0; jnode < sc.sizeOwned(); ++jnode) {
                gpfield(jnode, jvar) = rgpview_(f, jnode);
            }
            ++f;
        }
    }
    void unpack_3(const StructuredColumns& sc, Field& field, idx_t components) {
        auto gpfield = make_view<double, 3>(field);
        if (not components) {
            components = gpfield.shape(2);
        }
        for (idx_t jcomp = 0; jcomp < components; ++jcomp) {
            for (idx_t jlev = 0; jlev < gpfield.shape(1); ++jlev) {
                for (idx_t jnode = 0; jnode < sc.sizeOwned(); ++jnode) {
                    gpfield(jnode, jlev, jcomp) = rgpview_(f, jnode);
                }
                ++f;
            }
        }
    }
};

struct UnpackSpectral {
    const LocalView<double, 2>& rspecview_;
    size_t f;
    UnpackSpectral(const LocalView<double, 2>& rspecview): rspecview_(rspecview), f(0) {}

    void operator()(Field& field) {
        switch (field.rank()) {
            case 1:
                unpack_1(field);
                break;
            case 2:
                unpack_2(field);
                break;
            default:
                ATLAS_DEBUG_VAR(field.rank());
                ATLAS_NOTIMPLEMENTED;
                //break;
        }
    }

    void unpack_1(Field& field) {
        ArrayView<double, 1> spfield = make_view<double, 1>(field);

        for (idx_t jwave = 0; jwave < spfield.shape(0); ++jwave) {
            spfield(jwave) = rspecview_(jwave, f);
        }
        ++f;
    }
    void unpack_2(Field& field) {
        ArrayView<double, 2> spfield = make_view<double, 2>(field);

        const idx_t nvars = spfield.shape(1);

        for (idx_t jvar = 0; jvar < nvars; ++jvar) {
            for (idx_t jwave = 0; jwave < spfield.shape(0); ++jwave) {
                spfield(jwave, jvar) = rspecview_(jwave, f);
            }
            ++f;
        }
    }
};

}  // end anonymous namespace
}  // namespace atlas

namespace atlas {
namespace trans {

void TransIFS::assertCompatibleDistributions(const FunctionSpace& gp, const FunctionSpace& /*sp*/) const {
    std::string gp_dist = gp.distribution();
    if (gp_dist != "ectrans" &&  // distribution computed by TransPartitioner
        gp_dist != "serial" &&   // serial distribution always works
        gp_dist != "custom") {   // trust user that he knows what he is doing
        throw_Exception(gp.type() + " functionspace has unsupported distribution (" + gp_dist +
                            ") "
                            "to do spectral transforms. Please "
                            "partition grid with TransPartitioner",
                        Here());
    }
}

TransIFS::TransIFS(const Cache& cache, const Grid& grid, const long truncation, const eckit::Configuration& config):
    grid_(grid), cache_(cache.legendre().data()), cachesize_(cache.legendre().size()) {
    ctor(grid, truncation, config);
}

TransIFS::TransIFS(const Grid& grid, const long truncation, const eckit::Configuration& config):
    TransIFS(Cache(), grid, truncation, config) {}

TransIFS::TransIFS(const Grid& grid, const eckit::Configuration& config): TransIFS(grid, /*grid-only*/ -1, config) {}

int atlas::trans::TransIFS::handle() const {
    return trans_->handle;
}

int atlas::trans::TransIFS::ndgl() const {
    return trans_->ndgl;
}

int atlas::trans::TransIFS::nsmax() const {
    return trans_->nsmax;
}

int atlas::trans::TransIFS::ngptot() const {
    return trans_->ngptot;
}

int atlas::trans::TransIFS::ngptotg() const {
    return trans_->ngptotg;
}

int atlas::trans::TransIFS::ngptotmx() const {
    return trans_->ngptotmx;
}

int atlas::trans::TransIFS::nspec() const {
    return trans_->nspec;
}

int atlas::trans::TransIFS::nspec2() const {
    return trans_->nspec2;
}

int atlas::trans::TransIFS::nspec2g() const {
    return trans_->nspec2g;
}

int atlas::trans::TransIFS::nspec2mx() const {
    return trans_->nspec2mx;
}

int atlas::trans::TransIFS::n_regions_NS() const {
    return trans_->n_regions_NS;
}

int atlas::trans::TransIFS::n_regions_EW() const {
    return trans_->n_regions_EW;
}

int atlas::trans::TransIFS::nump() const {
    return trans_->nump;
}

int atlas::trans::TransIFS::nproc() const {
    return trans_->nproc;
}

int atlas::trans::TransIFS::myproc(int proc0) const {
    return trans_->myproc - 1 + proc0;
}

const int* atlas::trans::TransIFS::nloen(int& size) const {
    size = trans_->ndgl;
    ATLAS_ASSERT(trans_->nloen != nullptr);
    return trans_->nloen;
}

const int* atlas::trans::TransIFS::n_regions(int& size) const {
    size = trans_->n_regions_NS;
    ATLAS_ASSERT(trans_->n_regions != nullptr);
    return trans_->n_regions;
}

const int* atlas::trans::TransIFS::nfrstlat(int& size) const {
    size = trans_->n_regions_NS;
    if (trans_->nfrstlat == nullptr) {
        ::trans_inquire(trans_.get(), "nfrstlat");
    }
    return trans_->nfrstlat;
}

const int* atlas::trans::TransIFS::nlstlat(int& size) const {
    size = trans_->n_regions_NS;
    if (trans_->nlstlat == nullptr) {
        ::trans_inquire(trans_.get(), "nlstlat");
    }
    return trans_->nlstlat;
}

const int* atlas::trans::TransIFS::nptrfrstlat(int& size) const {
    size = trans_->n_regions_NS;
    if (trans_->nptrfrstlat == nullptr) {
        ::trans_inquire(trans_.get(), "nptrfrstlat");
    }
    return trans_->nptrfrstlat;
}

const int* atlas::trans::TransIFS::nsta(int& sizef2, int& sizef1) const {
    sizef1 = trans_->ndgl + trans_->n_regions_NS - 1;
    sizef2 = trans_->n_regions_EW;
    if (trans_->nsta == nullptr) {
        ::trans_inquire(trans_.get(), "nsta");
    }
    return trans_->nsta;
}

const int* atlas::trans::TransIFS::nonl(int& sizef2, int& sizef1) const {
    sizef1 = trans_->ndgl + trans_->n_regions_NS - 1;
    sizef2 = trans_->n_regions_EW;
    if (trans_->nonl == nullptr) {
        ::trans_inquire(trans_.get(), "nonl");
    }
    return trans_->nonl;
}

const int* atlas::trans::TransIFS::nmyms(int& size) const {
    size = trans_->nump;
    if (trans_->nmyms == nullptr) {
        ::trans_inquire(trans_.get(), "nmyms");
    }
    return trans_->nmyms;
}

const int* atlas::trans::TransIFS::nasm0(int& size) const {
    size = trans_->nsmax + 1;  // +1 because zeroth wave included
    if (trans_->nasm0 == nullptr) {
        ::trans_inquire(trans_.get(), "nasm0");
    }
    return trans_->nasm0;
}

const int* atlas::trans::TransIFS::nvalue(int& size) const {
    size = trans_->nspec2;
    if (trans_->nvalue == nullptr) {
        ::trans_inquire(trans_.get(), "nvalue");
    }
    return trans_->nvalue;
}

array::LocalView<int, 1> atlas::trans::TransIFS::nvalue() const {
    if (trans_->nvalue == nullptr) {
        ::trans_inquire(trans_.get(), "nvalue");
    }
    return array::LocalView<int, 1>(trans_->nvalue, array::make_shape(trans_->nspec2));
}

const int* atlas::trans::TransIFS::mvalue(int& size) const {
    size = trans_->nspec2;
    if( trans_->mvalue == nullptr ) {
        ::trans_inquire(trans_.get(), "mvalue");
    }
    return trans_->mvalue;
}

array::LocalView<int,1> atlas::trans::TransIFS::mvalue() const {
    if( trans_->mvalue == nullptr ) {
        ::trans_inquire(trans_.get(), "mvalue");
    }
    return array::LocalView<int,1> (trans_->mvalue, array::make_shape(trans_->nspec2));
}

array::LocalView<int, 1> atlas::trans::TransIFS::nasm0() const {
    if (trans_->nasm0 == nullptr) {
        ::trans_inquire(trans_.get(), "nasm0");
    }
    return array::LocalView<int, 1>(trans_->nasm0, array::make_shape(trans_->nsmax + 1));
}

array::LocalView<int, 1> atlas::trans::TransIFS::nmyms() const {
    if (trans_->nmyms == nullptr) {
        ::trans_inquire(trans_.get(), "nmyms");
    }
    return array::LocalView<int, 1>(trans_->nmyms, array::make_shape(trans_->nump));
}

array::LocalView<int, 2> atlas::trans::TransIFS::nonl() const {
    if (trans_->nonl == nullptr) {
        ::trans_inquire(trans_.get(), "nonl");
    }
    return array::LocalView<int, 2>(trans_->nonl,
                                    array::make_shape(trans_->n_regions_EW, trans_->ndgl + trans_->n_regions_NS - 1));
}

array::LocalView<int, 2> atlas::trans::TransIFS::nsta() const {
    if (trans_->nsta == nullptr) {
        ::trans_inquire(trans_.get(), "nsta");
    }
    return array::LocalView<int, 2>(trans_->nsta,
                                    array::make_shape(trans_->n_regions_EW, trans_->ndgl + trans_->n_regions_NS - 1));
}

array::LocalView<int, 1> atlas::trans::TransIFS::nptrfrstlat() const {
    if (trans_->nptrfrstlat == nullptr) {
        ::trans_inquire(trans_.get(), "nptrfrstlat");
    }
    return array::LocalView<int, 1>(trans_->nptrfrstlat, array::make_shape(trans_->n_regions_NS));
}

array::LocalView<int, 1> atlas::trans::TransIFS::nlstlat() const {
    if (trans_->nlstlat == nullptr) {
        ::trans_inquire(trans_.get(), "nlstlat");
    }
    return array::LocalView<int, 1>(trans_->nlstlat, array::make_shape(trans_->n_regions_NS));
}

array::LocalView<int, 1> atlas::trans::TransIFS::nfrstlat() const {
    if (trans_->nfrstlat == nullptr) {
        ::trans_inquire(trans_.get(), "nfrstlat");
    }
    return array::LocalView<int, 1>(trans_->nfrstlat, array::make_shape(trans_->n_regions_NS));
}

array::LocalView<int, 1> atlas::trans::TransIFS::n_regions() const {
    ATLAS_ASSERT(trans_->n_regions != nullptr);
    return array::LocalView<int, 1>(trans_->n_regions, array::make_shape(trans_->n_regions_NS));
}

array::LocalView<int, 1> atlas::trans::TransIFS::nloen() const {
    ATLAS_ASSERT(trans_->nloen != nullptr);
    return array::LocalView<int, 1>(trans_->nloen, array::make_shape(trans_->ndgl));
}


TransIFS::TransIFS(const Grid& grid, const Domain& domain, const long truncation, const eckit::Configuration& config):
    TransIFS(Cache(), grid, truncation, config) {
    ATLAS_ASSERT(domain.global());
}

TransIFS::TransIFS(const Cache& cache, const Grid& grid, const Domain& domain, const long truncation,
                   const eckit::Configuration& config):
    TransIFS(cache, grid, truncation, config) {
    ATLAS_ASSERT(domain.global());
}

TransIFS::TransIFS(const Grid& grid, const long truncation_x, const long truncation_y, const eckit::Configuration& config) :
    TransIFS(Cache(), grid, truncation_y, util::Config(config)("truncation_x",truncation_x)("truncation_y",truncation_y)) {
}

TransIFS::TransIFS(const Cache& cache, const Grid& grid, const Domain& domain, const long truncation_x, const long truncation_y,
                   const eckit::Configuration& config):
    TransIFS(grid, truncation_x, truncation_y, config) {
}


TransIFS::~TransIFS() = default;

int atlas::trans::TransIFS::truncation() const {
    return std::max(0, trans_->nsmax);
}

size_t TransIFS::nb_spectral_coefficients() const {
    return trans_->nspec2;
}

size_t TransIFS::nb_spectral_coefficients_global() const {
    return trans_->nspec2g;
}

const functionspace::Spectral& TransIFS::spectral() const {
    if (not spectral_) {
        spectral_ = functionspace::Spectral( truncation() );
    }
    return spectral_;
}

void TransIFS::ctor(const Grid& grid, long truncation, const eckit::Configuration& config) {
    trans_ = std::shared_ptr<::Trans_t>(new ::Trans_t, [](::Trans_t* p) {
        ::trans_delete(p);
        delete p;
    });
    if( auto rg = RegularGrid(grid); rg && not rg.domain().global() ){
        ctor_lam(grid, truncation, config);
        return;
    }
    ATLAS_ASSERT(grid.domain().global());
    ATLAS_ASSERT(not grid.projection());

    if (auto gg = GaussianGrid(grid)) {
        ctor_rgg(gg.ny(), gg.nx().data(), truncation, config);
        return;
    }
    if (auto ll = RegularLonLatGrid(grid)) {
        if (ll.standard() || ll.shifted()) {
            ctor_lonlat(ll.nx(), ll.ny(), truncation, config);
            return;
        }
    }
    throw_NotImplemented("Grid type not supported for Spectral Transforms", Here());
}

void TransIFS::ctor_rgg(const long nlat, const idx_t pl[], long truncation, const eckit::Configuration& config) {
    TransParameters p(*this, config);
    std::vector<int> nloen(nlat);
    for (long jlat = 0; jlat < nlat; ++jlat) {
        nloen[jlat] = pl[jlat];
    }
    TRANS_CHECK(::trans_new(trans_.get()));
    TRANS_CHECK(::trans_use_mpi(mpi::size() > 1));
    TRANS_CHECK(::trans_set_resol(trans_.get(), nlat, nloen.data()));
    if (truncation >= 0) {
        TRANS_CHECK(::trans_set_trunc(trans_.get(), truncation));
    }

    TRANS_CHECK(::trans_set_cache(trans_.get(), cache_, cachesize_));

    if (p.read_legendre().size() && mpi::size() == 1) {
        eckit::PathName file(p.read_legendre());
        if (not file.exists()) {
            std::stringstream msg;
            msg << "File " << file << " doesn't exist";
            throw_CantOpenFile(msg.str(), Here());
        }
        TRANS_CHECK(::trans_set_read(trans_.get(), file.asString().c_str()));
    }
    if (p.write_legendre().size() && mpi::size() == 1) {
        eckit::PathName file(p.write_legendre());
        TRANS_CHECK(::trans_set_write(trans_.get(), file.asString().c_str()));
    }

    trans_->fft    = p.fft();
    trans_->lsplit = p.split_y();
    trans_->flt    = p.flt();
    ATLAS_TRACE_SCOPE("trans_setup") { TRANS_CHECK(::trans_setup(trans_.get())); }
}

void TransIFS::ctor_lonlat(const long nlon, const long nlat, long truncation, const eckit::Configuration& config) {
    TransParameters p(*this, config);
    TRANS_CHECK(::trans_new(trans_.get()));
    TRANS_CHECK(::trans_use_mpi(mpi::size() > 1));
    TRANS_CHECK(::trans_set_resol_lonlat(trans_.get(), nlon, nlat));
    if (truncation >= 0) {
        TRANS_CHECK(::trans_set_trunc(trans_.get(), truncation));
    }
    TRANS_CHECK(::trans_set_cache(trans_.get(), cache_, cachesize_));

    if (p.read_legendre().size() && mpi::size() == 1) {
        eckit::PathName file(p.read_legendre());
        if (not file.exists()) {
            std::stringstream msg;
            msg << "File " << file << " doesn't exist";
            throw_CantOpenFile(msg.str(), Here());
        }
        TRANS_CHECK(::trans_set_read(trans_.get(), file.asString().c_str()));
    }
    if (p.write_legendre().size() && mpi::size() == 1) {
        eckit::PathName file(p.write_legendre());
        TRANS_CHECK(::trans_set_write(trans_.get(), file.asString().c_str()));
    }

    trans_->fft    = p.fft();
    trans_->lsplit = p.split_y();
    trans_->flt    = p.flt();

    TRANS_CHECK(::trans_setup(trans_.get()));
}

// --------------------------------------------------------------------------------------------

void TransIFS::ctor_lam(const RegularGrid& grid, long truncation, const eckit::Configuration& config) {
    ATLAS_ASSERT( grid.nx() > 0 );
    ATLAS_ASSERT( grid.ny() > 0 );

    TransParameters p(*this, config);

    util::Config partitioner_config;
    partitioner_config.set("type","checkerboard");
    partitioner_config.set("split_y", p.split_y());
    partitioner_config.set("split_x", true);
    auto partitioner = grid::Partitioner(partitioner_config);
    auto* checkerboard = dynamic_cast<grid::detail::partitioner::CheckerboardPartitioner*>(partitioner.get());
    int partitions_x = 1;
    int partitions_y = 1;
    if (checkerboard) {
        auto [px, py] = checkerboard->checkerboardDimensions(grid);
        partitions_x = px;
        partitions_y = py;
        TRANS_CHECK(::trans_set_nprgpew(partitions_x));
    }
    TRANS_CHECK(::trans_use_mpi(mpi::size() > 1));

    TRANS_CHECK(::trans_new(trans_.get()));
    const double dy = std::abs(grid.y(1)-grid.y(0));
    TRANS_CHECK(::trans_set_resol_lam(trans_.get(), grid.nx(), grid.ny(), grid.dx(), dy));
    trans_->lsplit = p.split_y();

    long truncation_x = -1;
    long truncation_y = -1;

    if (truncation >= 0) {
        if (config.has("truncation_x") && config.has("truncation_y")) {
            truncation_x = config.getInt("truncation_x");
            truncation_y = config.getInt("truncation_y");
        }
        else {
            truncation_y = truncation;
            truncation_x = grid.nx() * truncation / grid.ny();
        }
        TRANS_CHECK(::trans_set_trunc_lam(trans_.get(), truncation_x, truncation_y));

        trans_->fft = p.fft();
    }

    TRANS_CHECK(::trans_setup(trans_.get()));

    ATLAS_ASSERT(partitions_x == trans_->n_regions_EW);
    ATLAS_ASSERT(partitions_y == trans_->n_regions_NS);
    ATLAS_ASSERT(grid.size()  == trans_->ngptotg);
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans(const StructuredColumns& gp, const Field& gpfield, const Spectral& sp, Field& spfield,
                          const eckit::Configuration&) const {
    ATLAS_ASSERT(gpfield.functionspace() == 0 || functionspace::StructuredColumns(gpfield.functionspace()));
    ATLAS_ASSERT(spfield.functionspace() == 0 || functionspace::Spectral(spfield.functionspace()));

    assertCompatibleDistributions(gp, sp);

    if (compute_nfld(gpfield) != compute_nfld(spfield)) {
        throw_Exception("dirtrans: different number of gridpoint fields than spectral fields", Here());
    }
    if ((int)gpfield.shape(0) < ngptot()) {
        throw_Exception("dirtrans: slowest moving index must be >= ngptot", Here());
    }
    const int nfld = compute_nfld(gpfield);

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackStructuredColumns pack(rgpview);
        pack(gp, gpfield);
    }

    // Do transform
    {
        struct ::DirTrans_t transform = ::new_dirtrans(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();
        TRANS_CHECK(::trans_dirtrans(&transform));
    }

    // Unpack spectral
    {
        UnpackSpectral unpack(rspview);
        unpack(spfield);
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans(const functionspace::NodeColumns& gp, const Field& gpfield, const Spectral& sp,
                          Field& spfield, const eckit::Configuration& config) const {
    FieldSet gpfields;
    gpfields.add(gpfield);
    FieldSet spfields;
    spfields.add(spfield);
    __dirtrans(gp, gpfields, sp, spfields, config);
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans(const StructuredColumns& gp, const FieldSet& gpfields, const Spectral& sp, FieldSet& spfields,
                          const eckit::Configuration&) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const idx_t nfld = compute_nfld(gpfields);
    for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
        const Field& f = gpfields[jfld];
        ATLAS_ASSERT(f.functionspace() == 0 || functionspace::StructuredColumns(f.functionspace()));
    }

    const int trans_sp_nfld = compute_nfld(spfields);

    if (nfld != trans_sp_nfld) {
        throw_Exception("dirtrans: different number of gridpoint fields than spectral fields", Here());
    }
    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackStructuredColumns pack(rgpview);
        for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
            pack(gp, gpfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::DirTrans_t transform = ::new_dirtrans(trans_.get());
        transform.nscalar             = int(nfld);
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();

        TRANS_CHECK(::trans_dirtrans(&transform));
    }

    // Unpack the spectral fields
    {
        UnpackSpectral unpack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            unpack(spfields[jfld]);
        }
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans(const functionspace::NodeColumns& gp, const FieldSet& gpfields, const Spectral& sp,
                          FieldSet& spfields, const eckit::Configuration&) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const int nfld          = compute_nfld(gpfields);
    const int trans_sp_nfld = compute_nfld(spfields);

    if (nfld != trans_sp_nfld) {
        throw_Exception("dirtrans: different number of gridpoint fields than spectral fields", Here());
    }

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackNodeColumns pack(rgpview, gp);
        for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
            pack(gpfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::DirTrans_t transform = ::new_dirtrans(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();
        TRANS_CHECK(::trans_dirtrans(&transform));
    }

    // Unpack the spectral fields
    {
        UnpackSpectral unpack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            unpack(spfields[jfld]);
        }
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans_wind2vordiv(const functionspace::StructuredColumns& gp, const Field& gpwind,
                                      const Spectral& sp, Field& spvor, Field& spdiv,
                                      const eckit::Configuration&) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const size_t nfld = compute_nfld(spvor);
    if (spdiv.shape(0) != spvor.shape(0)) {
        throw_Exception("dirtrans:vorticity not compatible with divergence for 1st dimension", Here());
    }
    if (spdiv.shape(1) != spvor.shape(1)) {
        throw_Exception("dirtrans:vorticity not compatible with divergence for 2st dimension", Here());
    }
    const size_t nwindfld = compute_nfld(gpwind);
    if (nwindfld != 2 * nfld && nwindfld != 3 * nfld) {
        throw_Exception("dirtrans:wind field is not compatible with vorticity, divergence.", Here());
    }

    if (spdiv.shape(0) != nspec2()) {
        std::stringstream msg;
        msg << "dirtrans: Spectral vorticity and divergence have wrong dimension: "
               "nspec2 "
            << spdiv.shape(0) << " should be " << nspec2();
        throw_Exception(msg.str(), Here());
    }

    if (spvor.size() == 0) {
        throw_Exception("dirtrans: spectral vorticity field is empty.");
    }
    if (spdiv.size() == 0) {
        throw_Exception("dirtrans: spectral divergence field is empty.");
    }

    // Arrays Trans expects
    std::vector<double> rgp(2 * nfld * ngptot());
    std::vector<double> rspvor(nspec2() * nfld);
    std::vector<double> rspdiv(nspec2() * nfld);
    auto rgpview    = LocalView<double, 2>(rgp.data(),    make_shape(2 * nfld, ngptot()));
    auto rspvorview = LocalView<double, 2>(rspvor.data(), make_shape(nspec2(), nfld));
    auto rspdivview = LocalView<double, 2>(rspdiv.data(), make_shape(nspec2(), nfld));

    std::vector<double> rmeanu;
    std::vector<double> rmeanv;

    // Pack gridpoints
    {
        PackStructuredColumns pack(rgpview);
        int wind_components = 2;
        pack(gpwind.functionspace(), gpwind, wind_components);
    }

    // Do transform
    {
        struct ::DirTrans_t transform = ::new_dirtrans(trans_.get());
        transform.nvordiv             = int(nfld);
        transform.rgp                 = rgp.data();
        transform.rspvor              = rspvor.data();
        transform.rspdiv              = rspdiv.data();
        if( trans_->llam ) {
            rmeanu.resize(nfld);
            rmeanv.resize(nfld);
            transform.rmeanu              = rmeanu.data();
            transform.rmeanv              = rmeanv.data();
        }
        ATLAS_ASSERT(transform.rspvor);
        ATLAS_ASSERT(transform.rspdiv);
        TRANS_CHECK(::trans_dirtrans(&transform));
    }

    // Pack spectral fields
    UnpackSpectral unpack_vor(rspvorview);
    unpack_vor(spvor);
    UnpackSpectral unpack_div(rspdivview);
    unpack_div(spdiv);

    ATLAS_DEBUG_VAR(rmeanu);
    ATLAS_DEBUG_VAR(rmeanv);
    if( trans_->llam ) {
        spvor.metadata().set("rmeanu",rmeanu);
        spvor.metadata().set("rmeanv",rmeanv);
    }
}

// -----------------------------------------------------------------------------------------------

void TransIFS::__dirtrans_wind2vordiv(const functionspace::NodeColumns& gp, const Field& gpwind, const Spectral& sp,
                                      Field& spvor, Field& spdiv, const eckit::Configuration&) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const size_t nfld = compute_nfld(spvor);
    if (spdiv.shape(0) != spvor.shape(0)) {
        throw_Exception("dirtrans: vorticity not compatible with divergence for 1st dimension.", Here());
    }
    if (spdiv.shape(1) != spvor.shape(1)) {
        throw_Exception("dirtrans: vorticity not compatible with divergence for 2nd dimension.", Here());
    }
    const size_t nwindfld = compute_nfld(gpwind);
    if (nwindfld != 2 * nfld && nwindfld != 3 * nfld) {
        throw_Exception("dirtrans: wind field is not compatible with vorticity, divergence.", Here());
    }

    if (spdiv.shape(0) != nspec2()) {
        std::stringstream msg;
        msg << "dirtrans: Spectral vorticity and divergence have wrong dimension: "
               "nspec2 "
            << spdiv.shape(0) << " should be " << nspec2();
        throw_Exception(msg.str(), Here());
    }

    if (spvor.size() == 0) {
        throw_Exception("dirtrans: spectral vorticity field is empty.");
    }
    if (spdiv.size() == 0) {
        throw_Exception("dirtrans: spectral divergence field is empty.");
    }

    // Arrays Trans expects
    std::vector<double> rgp(2 * nfld * ngptot());
    std::vector<double> rspvor(nspec2() * nfld);
    std::vector<double> rspdiv(nspec2() * nfld);
    auto rgpview    = LocalView<double, 2>(rgp.data(), make_shape(2 * nfld, ngptot()));
    auto rspvorview = LocalView<double, 2>(rspvor.data(), make_shape(nspec2(), nfld));
    auto rspdivview = LocalView<double, 2>(rspdiv.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackNodeColumns pack(rgpview, gp);
        int wind_components = 2;
        pack(gpwind, wind_components);
    }

    // Do transform
    {
        struct ::DirTrans_t transform = ::new_dirtrans(trans_.get());
        transform.nvordiv             = int(nfld);
        transform.rgp                 = rgp.data();
        transform.rspvor              = rspvor.data();
        transform.rspdiv              = rspdiv.data();

        ATLAS_ASSERT(transform.rspvor);
        ATLAS_ASSERT(transform.rspdiv);
        TRANS_CHECK(::trans_dirtrans(&transform));
    }

    // Pack spectral fields
    UnpackSpectral unpack_vor(rspvorview);
    unpack_vor(spvor);
    UnpackSpectral unpack_div(rspdivview);
    unpack_div(spdiv);
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans_adj(const Spectral& sp, const FieldSet& spfields,
                              const StructuredColumns& gp, FieldSet& gpfields,
                              const eckit::Configuration&) const {
#if ATLAS_HAVE_ECTRANS
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const idx_t nfld = compute_nfld(gpfields);
    for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
        const Field& f = gpfields[jfld];
        ATLAS_ASSERT(f.functionspace() == 0 || functionspace::StructuredColumns(f.functionspace()));
    }

    const int trans_sp_nfld = compute_nfld(spfields);

    if (nfld != trans_sp_nfld) {
        throw_Exception("dirtrans_adj: different number of gridpoint fields than spectral fields", Here());
    }
    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            pack(spfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::DirTransAdj_t transform = ::new_dirtrans_adj(trans_.get());
        transform.nscalar             = int(nfld);
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();

        TRANS_CHECK(::trans_dirtrans_adj(&transform));
    }

    //  Unpack the gridpoint fields
    {
        UnpackStructuredColumns unpack(rgpview);
        for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
            unpack(gp, gpfields[jfld]);
        }
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans_adj( const Spectral& sp, const Field& spfield,
                               const NodeColumns& gp, Field& gpfield,
                               const eckit::Configuration& ) const {
#if ATLAS_HAVE_ECTRANS
    ATLAS_ASSERT(gpfield.functionspace() == 0 || functionspace::NodeColumns(gpfield.functionspace()));
    ATLAS_ASSERT(spfield.functionspace() == 0 || functionspace::Spectral(spfield.functionspace()));

    assertCompatibleDistributions(gp, sp);

    if (compute_nfld(gpfield) != compute_nfld(spfield)) {
        throw_Exception("dirtrans: different number of gridpoint fields than spectral fields", Here());
    }
    if ((int)gpfield.shape(0) < ngptot()) {
        throw_Exception("dirtrans: slowest moving index must be >= ngptot", Here());
    }
    const int nfld = compute_nfld(gpfield);

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        pack(spfield);
    }

    // Do transform
    {
        struct ::DirTransAdj_t transform = ::new_dirtrans_adj(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();
        TRANS_CHECK(::trans_dirtrans_adj(&transform));
    }

    // Unpack gridpoint fields
    {
        UnpackNodeColumns unpack(rgpview, gp);
        unpack(gpfield);
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans_adj( const Spectral& sp, const Field& spfield,
                               const StructuredColumns& gp, Field& gpfield,
                               const eckit::Configuration& ) const {
#if ATLAS_HAVE_ECTRANS
    ATLAS_ASSERT(gpfield.functionspace() == 0 || functionspace::StructuredColumns(gpfield.functionspace()));
    ATLAS_ASSERT(spfield.functionspace() == 0 || functionspace::Spectral(spfield.functionspace()));

    assertCompatibleDistributions(gp, sp);

    if (compute_nfld(gpfield) != compute_nfld(spfield)) {
        throw_Exception("dirtrans: different number of gridpoint fields than spectral fields", Here());
    }
    if ((int)gpfield.shape(0) < ngptot()) {
        throw_Exception("dirtrans: slowest moving index must be >= ngptot", Here());
    }
    const int nfld = compute_nfld(gpfield);

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        pack(spfield);
    }

    // Do transform
    {
        struct ::DirTransAdj_t transform = ::new_dirtrans_adj(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();
        TRANS_CHECK(::trans_dirtrans_adj(&transform));
    }

    // Unpack gridpoint fields
    {
        UnpackStructuredColumns unpack(rgpview);
        unpack(gp, gpfield);
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}


// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans_adj(const Spectral& sp, const FieldSet& spfields,
                              const NodeColumns& gp, FieldSet& gpfields,
                              const eckit::Configuration&) const {
#if ATLAS_HAVE_ECTRANS
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const idx_t nfld = compute_nfld(gpfields);
    for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
        const Field& f = gpfields[jfld];
        ATLAS_ASSERT(f.functionspace() == 0 || functionspace::NodeColumns(f.functionspace()));
    }

    const int trans_sp_nfld = compute_nfld(spfields);

    if (nfld != trans_sp_nfld) {
        throw_Exception("dirtrans_adj: different number of gridpoint fields than spectral fields", Here());
    }
    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            pack(spfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::DirTransAdj_t transform = ::new_dirtrans_adj(trans_.get());
        transform.nscalar             = int(nfld);
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();

        TRANS_CHECK(::trans_dirtrans_adj(&transform));
    }

    //  Unpack the gridpoint fields
    {
        UnpackNodeColumns unpack(rgpview, gp);
        for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
            unpack(gpfields[jfld]);
        }
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans_wind2vordiv_adj(const Spectral& sp, const Field& spvor, const Field& spdiv,
                                          const functionspace::StructuredColumns& gp, Field& gpwind,
                                          const eckit::Configuration&) const {
#if ATLAS_HAVE_ECTRANS
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const size_t nfld = compute_nfld(spvor);
    if (spdiv.shape(0) != spvor.shape(0)) {
        throw_Exception("dirtrans_adj:vorticity not compatible with divergence for 1st dimension", Here());
    }
    if (spdiv.shape(1) != spvor.shape(1)) {
        throw_Exception("dirtrans_adj:vorticity not compatible with divergence for 2st dimension", Here());
    }
    const size_t nwindfld = compute_nfld(gpwind);
    if (nwindfld != 2 * nfld && nwindfld != 3 * nfld) {
        throw_Exception("dirtrans_adj:wind field is not compatible with vorticity, divergence.", Here());
    }

    if (spdiv.shape(0) != nspec2()) {
        std::stringstream msg;
        msg << "dirtrans_adj: Spectral vorticity and divergence have wrong dimension: "
               "nspec2 "
            << spdiv.shape(0) << " should be " << nspec2();
        throw_Exception(msg.str(), Here());
    }

    if (spvor.size() == 0) {
        throw_Exception("dirtrans_adj: spectral vorticity field is empty.");
    }
    if (spdiv.size() == 0) {
        throw_Exception("dirtrans_adj: spectral divergence field is empty.");
    }

    // Arrays Trans expects
    std::vector<double> rgp(2 * nfld * ngptot());
    std::vector<double> rspvor(nspec2() * nfld);
    std::vector<double> rspdiv(nspec2() * nfld);
    auto rgpview    = LocalView<double, 2>(rgp.data(), make_shape(2 * nfld, ngptot()));
    auto rspvorview = LocalView<double, 2>(rspvor.data(), make_shape(nspec2(), nfld));
    auto rspdivview = LocalView<double, 2>(rspdiv.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    PackSpectral pack_vor(rspvorview);
    pack_vor(spvor);
    PackSpectral pack_div(rspdivview);
    pack_div(spdiv);

    // Do transform
    {
        struct ::DirTransAdj_t transform = ::new_dirtrans_adj(trans_.get());
        transform.nvordiv             = int(nfld);
        transform.rgp                 = rgp.data();
        transform.rspvor              = rspvor.data();
        transform.rspdiv              = rspdiv.data();

        ATLAS_ASSERT(transform.rspvor);
        ATLAS_ASSERT(transform.rspdiv);
        TRANS_CHECK(::trans_dirtrans_adj(&transform));
    }

    // Unpack the gridpoint fields
    {
        UnpackStructuredColumns unpack(rgpview);
        int wind_components = 2;
        unpack(gp, gpwind, wind_components);
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

// --------------------------------------------------------------------------------------------

void TransIFS::__dirtrans_wind2vordiv_adj(const Spectral& sp, const Field& spvor, const Field& spdiv,
                                          const functionspace::NodeColumns& gp, Field& gpwind,
                                          const eckit::Configuration&) const {
#if ATLAS_HAVE_ECTRANS
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const size_t nfld = compute_nfld(spvor);
    if (spdiv.shape(0) != spvor.shape(0)) {
        throw_Exception("dirtrans_adj:vorticity not compatible with divergence for 1st dimension", Here());
    }
    if (spdiv.shape(1) != spvor.shape(1)) {
        throw_Exception("dirtrans_adj:vorticity not compatible with divergence for 2st dimension", Here());
    }
    const size_t nwindfld = compute_nfld(gpwind);
    if (nwindfld != 2 * nfld && nwindfld != 3 * nfld) {
        throw_Exception("dirtrans_adj:wind field is not compatible with vorticity, divergence.", Here());
    }

    if (spdiv.shape(0) != nspec2()) {
        std::stringstream msg;
        msg << "dirtrans_adj: Spectral vorticity and divergence have wrong dimension: "
               "nspec2 "
            << spdiv.shape(0) << " should be " << nspec2();
        throw_Exception(msg.str(), Here());
    }

    if (spvor.size() == 0) {
        throw_Exception("dirtrans_adj: spectral vorticity field is empty.");
    }
    if (spdiv.size() == 0) {
        throw_Exception("dirtrans_adj: spectral divergence field is empty.");
    }

    // Arrays Trans expects
    std::vector<double> rgp(2 * nfld * ngptot());
    std::vector<double> rspvor(nspec2() * nfld);
    std::vector<double> rspdiv(nspec2() * nfld);
    auto rgpview    = LocalView<double, 2>(rgp.data(), make_shape(2 * nfld, ngptot()));
    auto rspvorview = LocalView<double, 2>(rspvor.data(), make_shape(nspec2(), nfld));
    auto rspdivview = LocalView<double, 2>(rspdiv.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    PackSpectral pack_vor(rspvorview);
    pack_vor(spvor);
    PackSpectral pack_div(rspdivview);
    pack_div(spdiv);

    // Do transform
    {
        struct ::DirTransAdj_t transform = ::new_dirtrans_adj(trans_.get());
        transform.nvordiv             = int(nfld);
        transform.rgp                 = rgp.data();
        transform.rspvor              = rspvor.data();
        transform.rspdiv              = rspdiv.data();

        ATLAS_ASSERT(transform.rspvor);
        ATLAS_ASSERT(transform.rspdiv);
        TRANS_CHECK(::trans_dirtrans_adj(&transform));
    }

    // Unpack the gridpoint fields
    {
        UnpackNodeColumns unpack(rgpview, gp);
        int wind_components = 2;
        unpack(gpwind, wind_components);
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}



// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans(const functionspace::Spectral& sp, const Field& spfield,
                          const functionspace::StructuredColumns& gp, Field& gpfield,
                          const eckit::Configuration& config) const {
    assertCompatibleDistributions(gp, sp);

    ATLAS_ASSERT(gpfield.functionspace() == 0 || functionspace::StructuredColumns(gpfield.functionspace()));
    ATLAS_ASSERT(spfield.functionspace() == 0 || functionspace::Spectral(spfield.functionspace()));
    if (compute_nfld(gpfield) != compute_nfld(spfield)) {
        throw_Exception("invtrans: different number of gridpoint fields than spectral fields", Here());
    }
    if ((int)gpfield.shape(0) < ngptot()) {
        throw_Exception("invtrans: slowest moving index must be >= ngptot", Here());
    }
    const int nfld = compute_nfld(gpfield);

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        pack(spfield);
    }

    // Do transform
    {
        struct ::InvTrans_t transform = ::new_invtrans(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();
        TRANS_CHECK(::trans_invtrans(&transform));
    }

    // Unpack gridpoint fields
    {
        UnpackStructuredColumns unpack(rgpview);
        unpack(gp, gpfield);
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans(const Spectral& sp, const Field& spfield, const functionspace::NodeColumns& gp,
                          Field& gpfield, const eckit::Configuration& config) const {
    FieldSet spfields;
    spfields.add(spfield);
    FieldSet gpfields;
    gpfields.add(gpfield);
    __invtrans(sp, spfields, gp, gpfields, config);
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans(const functionspace::Spectral& sp, const FieldSet& spfields,
                          const functionspace::StructuredColumns& gp, FieldSet& gpfields,
                          const eckit::Configuration& config) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const int nfld = compute_nfld(gpfields);
    for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
        const Field& f = gpfields[jfld];
        ATLAS_ASSERT(f.functionspace() == 0 || functionspace::StructuredColumns(f.functionspace()));
    }

    const int nb_spectral_fields = compute_nfld(spfields);

    if (nfld != nb_spectral_fields) {
        std::stringstream msg;
        msg << "invtrans: different number of gridpoint fields than spectral fields"
            << "[ " << nfld << " != " << nb_spectral_fields << " ]";
        throw_Exception(msg.str(), Here());
    }

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            pack(spfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::InvTrans_t transform = ::new_invtrans(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();

        TRANS_CHECK(::trans_invtrans(&transform));
    }

    // Unpack the gridpoint fields
    {
        UnpackStructuredColumns unpack(rgpview);
        for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
            unpack(gp, gpfields[jfld]);
        }
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans(const Spectral& sp, const FieldSet& spfields, const functionspace::NodeColumns& gp,
                          FieldSet& gpfields, const eckit::Configuration& config) const {
    assertCompatibleDistributions(gp, sp);


    // Count total number of fields and do sanity checks
    const int nfld               = compute_nfld(gpfields);
    const int nb_spectral_fields = compute_nfld(spfields);

    if (nfld != nb_spectral_fields) {
        throw_Exception("invtrans: different number of gridpoint fields than spectral fields", Here());
    }

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            pack(spfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::InvTrans_t transform = ::new_invtrans(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();
        TRANS_CHECK(::trans_invtrans(&transform));
    }

    // Unpack the gridpoint fields
    {
        UnpackNodeColumns unpack(rgpview, gp);
        for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
            unpack(gpfields[jfld]);
        }
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_grad(const Spectral& sp, const Field& spfield, const functionspace::StructuredColumns& gp,
                               Field& gradfield, const eckit::Configuration& config) const {
    FieldSet spfields;
    spfields.add(spfield);
    FieldSet gradfields;
    gradfields.add(gradfield);
    __invtrans_grad(sp, spfields, gp, gradfields, config);
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_grad(const Spectral& sp, const Field& spfield, const functionspace::NodeColumns& gp,
                               Field& gradfield, const eckit::Configuration& config) const {
    FieldSet spfields;
    spfields.add(spfield);
    FieldSet gradfields;
    gradfields.add(gradfield);
    __invtrans_grad(sp, spfields, gp, gradfields, config);
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_grad(const Spectral& sp, const FieldSet& spfields, const functionspace::StructuredColumns& gp,
                               FieldSet& gradfields, const eckit::Configuration& config) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const int nb_gridpoint_field = compute_nfld(gradfields);
    const int nfld               = compute_nfld(spfields);

    if (nb_gridpoint_field != 2 * nfld) {  // factor 2 because N-S and E-W derivatives
        throw_Exception(
            "invtrans_grad: different number of gridpoint "
            "fields than spectral fields",
            Here());
    }

    // Arrays Trans expects
    // Allocate space for
    std::vector<double> rgp(3 * nfld * ngptot());  // (scalars) + (NS ders) + (EW ders)
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(3 * nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));


    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            pack(spfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::InvTrans_t transform = ::new_invtrans(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();
        transform.lscalarders         = true;

        TRANS_CHECK(::trans_invtrans(&transform));
    }

    // Unpack the gridpoint fields
    {
        int f = nfld;  // skip to where derivatives start
        for (idx_t dim = 0; dim < 2; ++dim) {
            for (idx_t jfld = 0; jfld < gradfields.size(); ++jfld) {
                const idx_t nb_nodes = StructuredColumns(gradfields[jfld].functionspace()).sizeOwned();
                const idx_t nlev     = gradfields[jfld].levels();
                if (nlev) {
                    auto field = make_view<double, 3>(gradfields[jfld]);
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        for (idx_t jnode = 0; jnode < nb_nodes; ++jnode) {
                            field(jnode, jlev, 1 - dim) = rgpview(f, jnode);
                        }
                    }
                }
                else {
                    auto field = make_view<double, 2>(gradfields[jfld]);
                    for (idx_t jnode = 0; jnode < nb_nodes; ++jnode) {
                        field(jnode, 1 - dim) = rgpview(f, jnode);
                    }
                }
                ++f;
            }
        }
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_grad(const Spectral& sp, const FieldSet& spfields, const functionspace::NodeColumns& gp,
                               FieldSet& gradfields, const eckit::Configuration& config) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const int nb_gridpoint_field = compute_nfld(gradfields);
    const int nfld               = compute_nfld(spfields);

    if (nb_gridpoint_field != 2 * nfld) {  // factor 2 because N-S and E-W derivatives
        throw_Exception(
            "invtrans_grad: different number of gridpoint "
            "fields than spectral fields",
            Here());
    }

    // Arrays Trans expects
    // Allocate space for
    std::vector<double> rgp(3 * nfld * ngptot());  // (scalars) + (NS ders) + (EW ders)
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(3 * nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));


    // Pack spectral fields
    {
        PackSpectral pack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            pack(spfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::InvTrans_t transform = ::new_invtrans(trans_.get());
        transform.nscalar             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspscalar           = rsp.data();
        transform.lscalarders         = true;

        TRANS_CHECK(::trans_invtrans(&transform));
    }

    // Unpack the gridpoint fields
    {
        mesh::IsGhostNode is_ghost(gp.nodes());
        int f = nfld;  // skip to where derivatives start
        for (idx_t dim = 0; dim < 2; ++dim) {
            for (idx_t jfld = 0; jfld < gradfields.size(); ++jfld) {
                const idx_t nb_nodes = gradfields[jfld].shape(0);
                const idx_t nlev     = gradfields[jfld].levels();
                if (nlev) {
                    auto field = make_view<double, 3>(gradfields[jfld]);
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        int n = 0;
                        for (idx_t jnode = 0; jnode < nb_nodes; ++jnode) {
                            if (!is_ghost(jnode)) {
                                field(jnode, jlev, 1 - dim) = rgpview(f, n);
                                ++n;
                            }
                        }
                        ATLAS_ASSERT(n == ngptot());
                    }
                }
                else {
                    auto field = make_view<double, 2>(gradfields[jfld]);
                    int n      = 0;
                    for (idx_t jnode = 0; jnode < nb_nodes; ++jnode) {
                        if (!is_ghost(jnode)) {
                            field(jnode, 1 - dim) = rgpview(f, n);
                            ++n;
                        }
                    }
                    ATLAS_ASSERT(n == ngptot());
                }
                ++f;
            }
        }
    }
}

// --------------------------------------------------------------------------------------------
void TransIFS::__invtrans_vordiv2wind(const Spectral& sp, const Field& spvor, const Field& spdiv,
                                      const functionspace::StructuredColumns& gp, Field& gpwind,
                                      const eckit::Configuration&) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const int nfld = compute_nfld(spvor);
    if (spdiv.shape(0) != spvor.shape(0)) {
        throw_Exception("invtrans: vorticity not compatible with divergence.", Here());
    }
    if (spdiv.shape(1) != spvor.shape(1)) {
        throw_Exception("invtrans: vorticity not compatible with divergence.", Here());
    }
    const int nwindfld = compute_nfld(gpwind);
    if (nwindfld != 2 * nfld && nwindfld != 3 * nfld) {
        throw_Exception("invtrans: wind field is not compatible with vorticity, divergence.", Here());
    }

    if (spdiv.shape(0) != nspec2()) {
        std::stringstream msg;
        msg << "invtrans: Spectral vorticity and divergence have wrong dimension: "
               "nspec2 "
            << spdiv.shape(0) << " should be " << nspec2();
        throw_Exception(msg.str(), Here());
    }

    ATLAS_ASSERT(spvor.rank() == 2);
    ATLAS_ASSERT(spdiv.rank() == 2);
    if (spvor.size() == 0) {
        throw_Exception("invtrans: spectral vorticity field is empty.");
    }
    if (spdiv.size() == 0) {
        throw_Exception("invtrans: spectral divergence field is empty.");
    }

    // Arrays Trans expects
    std::vector<double> rgp(2 * nfld * ngptot());
    std::vector<double> rspvor(nspec2() * nfld);
    std::vector<double> rspdiv(nspec2() * nfld);
    auto rgpview    = LocalView<double, 2>(rgp.data(), make_shape(2 * nfld, ngptot()));
    auto rspvorview = LocalView<double, 2>(rspvor.data(), make_shape(nspec2(), nfld));
    auto rspdivview = LocalView<double, 2>(rspdiv.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    PackSpectral pack_vor(rspvorview);
    pack_vor(spvor);
    PackSpectral pack_div(rspdivview);
    pack_div(spdiv);

    // Do transform
    {
        std::vector<double> rmeanu;
        std::vector<double> rmeanv;

        struct ::InvTrans_t transform = ::new_invtrans(trans_.get());
        transform.nvordiv             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspvor              = rspvor.data();
        transform.rspdiv              = rspdiv.data();
        if (trans_->llam) {
            rmeanu = spvor.metadata().getDoubleVector("rmeanu");
            rmeanv = spvor.metadata().getDoubleVector("rmeanv");
            transform.rmeanu = rmeanu.data();
            transform.rmeanv = rmeanv.data();
        }

        ATLAS_ASSERT(transform.rspvor);
        ATLAS_ASSERT(transform.rspdiv);
        TRANS_CHECK(::trans_invtrans(&transform));
    }

    // Unpack the gridpoint fields
    {
        UnpackStructuredColumns unpack(rgpview);
        int wind_components = 2;
        unpack(gp, gpwind, wind_components);
    }
}

//---------------------------------------------------------------------------------------------

void TransIFS::__invtrans_vordiv2wind(const Spectral& sp, const Field& spvor, const Field& spdiv,
                                      const functionspace::NodeColumns& gp, Field& gpwind,
                                      const eckit::Configuration&) const {
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const int nfld = compute_nfld(spvor);
    if (spdiv.shape(0) != spvor.shape(0)) {
        throw_Exception("invtrans: vorticity not compatible with divergence.", Here());
    }
    if (spdiv.shape(1) != spvor.shape(1)) {
        throw_Exception("invtrans: vorticity not compatible with divergence.", Here());
    }
    const int nwindfld = compute_nfld(gpwind);
    if (nwindfld != 2 * nfld && nwindfld != 3 * nfld) {
        throw_Exception("invtrans: wind field is not compatible with vorticity, divergence.", Here());
    }

    if (spdiv.shape(0) != nspec2()) {
        std::stringstream msg;
        msg << "invtrans: Spectral vorticity and divergence have wrong dimension: "
               "nspec2 "
            << spdiv.shape(0) << " should be " << nspec2();
        throw_Exception(msg.str(), Here());
    }

    ATLAS_ASSERT(spvor.rank() == 2);
    ATLAS_ASSERT(spdiv.rank() == 2);
    if (spvor.size() == 0) {
        throw_Exception("invtrans: spectral vorticity field is empty.");
    }
    if (spdiv.size() == 0) {
        throw_Exception("invtrans: spectral divergence field is empty.");
    }

    // Arrays Trans expects
    std::vector<double> rgp(2 * nfld * ngptot());
    std::vector<double> rspvor(nspec2() * nfld);
    std::vector<double> rspdiv(nspec2() * nfld);
    auto rgpview    = LocalView<double, 2>(rgp.data(), make_shape(2 * nfld, ngptot()));
    auto rspvorview = LocalView<double, 2>(rspvor.data(), make_shape(nspec2(), nfld));
    auto rspdivview = LocalView<double, 2>(rspdiv.data(), make_shape(nspec2(), nfld));

    // Pack spectral fields
    PackSpectral pack_vor(rspvorview);
    pack_vor(spvor);
    PackSpectral pack_div(rspdivview);
    pack_div(spdiv);

    // Do transform
    {
        struct ::InvTrans_t transform = ::new_invtrans(trans_.get());
        transform.nvordiv             = nfld;
        transform.rgp                 = rgp.data();
        transform.rspvor              = rspvor.data();
        transform.rspdiv              = rspdiv.data();

        ATLAS_ASSERT(transform.rspvor);
        ATLAS_ASSERT(transform.rspdiv);
        TRANS_CHECK(::trans_invtrans(&transform));
    }

    // Unpack the gridpoint fields
    {
        UnpackNodeColumns unpack(rgpview, gp);
        int wind_components = 2;
        unpack(gpwind, wind_components);
    }
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_adj(const functionspace::Spectral& sp, Field& spfield,
                              const functionspace::StructuredColumns& gp, const Field& gpfield,
                              const eckit::Configuration& config) const {
#if ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ)

    ATLAS_ASSERT(gpfield.functionspace() == 0 || functionspace::StructuredColumns(gpfield.functionspace()));
    ATLAS_ASSERT(spfield.functionspace() == 0 || functionspace::Spectral(spfield.functionspace()));

    assertCompatibleDistributions(gp, sp);

    if (compute_nfld(gpfield) != compute_nfld(spfield)) {
        throw_Exception("dirtrans: different number of gridpoint fields than spectral fields", Here());
    }
    if ((int)gpfield.shape(0) < ngptot()) {
        throw_Exception("dirtrans: slowest moving index must be >= ngptot", Here());
    }
    const int nfld = compute_nfld(gpfield);

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackStructuredColumns pack(rgpview);
        pack(gp, gpfield);
    }

    // Do transform
    {
        struct ::InvTransAdj_t transform = ::new_invtrans_adj(trans_.get());
        transform.nscalar                = nfld;
        transform.rgp                    = rgp.data();
        transform.rspscalar              = rsp.data();
        TRANS_CHECK(::trans_invtrans_adj(&transform));
    }

    // Unpack spectral
    {
        UnpackSpectral unpack(rspview);
        unpack(spfield);
    }

#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_adj(const Spectral& sp, Field& spfield, const functionspace::NodeColumns& gp,
                              const Field& gpfield, const eckit::Configuration& config) const {
    FieldSet spfields;
    spfields.add(spfield);
    FieldSet gpfields;
    gpfields.add(gpfield);
    __invtrans_adj(sp, spfields, gp, gpfields, config);
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_adj(const functionspace::Spectral& sp, FieldSet& spfields,
                              const functionspace::StructuredColumns& gp, const FieldSet& gpfields,
                              const eckit::Configuration& config) const {
#if ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ)

    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const idx_t nfld = compute_nfld(gpfields);
    for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
        const Field& f = gpfields[jfld];
        ATLAS_ASSERT(f.functionspace() == 0 || functionspace::StructuredColumns(f.functionspace()));
    }

    const int trans_sp_nfld = compute_nfld(spfields);

    if (nfld != trans_sp_nfld) {
        throw_Exception("invtrans_adj: different number of gridpoint fields than spectral fields", Here());
    }

    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackStructuredColumns pack(rgpview);
        for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
            pack(gp, gpfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::InvTransAdj_t transform = ::new_invtrans_adj(trans_.get());
        transform.nscalar                = int(nfld);
        transform.rgp                    = rgp.data();
        transform.rspscalar              = rsp.data();

        TRANS_CHECK(::trans_invtrans_adj(&transform));
    }

    // Unpack the spectral fields
    {
        UnpackSpectral unpack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            unpack(spfields[jfld]);
        }
    }

#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_adj(const Spectral& sp, FieldSet& spfields, const functionspace::NodeColumns& gp,
                              const FieldSet& gpfields, const eckit::Configuration& config) const {
#if ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ)

    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const idx_t nfld = compute_nfld(gpfields);
    for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
        const Field& f = gpfields[jfld];
        ATLAS_ASSERT(f.functionspace() == 0 || functionspace::StructuredColumns(f.functionspace()));
    }

    const int trans_sp_nfld = compute_nfld(spfields);

    if (nfld != trans_sp_nfld) {
        throw_Exception("invtrans_adj: different number of gridpoint fields than spectral fields", Here());
    }
    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackNodeColumns pack(rgpview, gp);
        for (idx_t jfld = 0; jfld < gpfields.size(); ++jfld) {
            pack(gpfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::InvTransAdj_t transform = ::new_invtrans_adj(trans_.get());
        transform.nscalar                = int(nfld);
        transform.rgp                    = rgp.data();
        transform.rspscalar              = rsp.data();

        TRANS_CHECK(::trans_invtrans_adj(&transform));
    }

    // Unpack the spectral fields
    {
        UnpackSpectral unpack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            unpack(spfields[jfld]);
        }
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

//---------------------------------------------------------------------------------------------

void TransIFS::__invtrans_grad_adj(const Spectral& sp, Field& spfield, const functionspace::StructuredColumns& gp,
                                   const Field& gradfield, const eckit::Configuration& config) const {
    FieldSet spfields;
    spfields.add(spfield);
    FieldSet gradfields;
    gradfields.add(gradfield);
    __invtrans_grad_adj(sp, spfields, gp, gradfields, config);
}

//---------------------------------------------------------------------------------------------

void TransIFS::__invtrans_grad_adj(const Spectral& sp, Field& spfield, const functionspace::NodeColumns& gp,
                                   const Field& gradfield, const eckit::Configuration& config) const {
    FieldSet spfields;
    spfields.add(spfield);
    FieldSet gradfields;
    gradfields.add(gradfield);
    __invtrans_grad_adj(sp, spfields, gp, gradfields, config);
}

//---------------------------------------------------------------------------------------------

void TransIFS::__invtrans_grad_adj(const Spectral& sp, FieldSet& spfields, const functionspace::StructuredColumns& gp,
                                   const FieldSet& gradfields, const eckit::Configuration& config) const {
#if ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ)
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const idx_t nfld = compute_nfld(gradfields);
    for (idx_t jfld = 0; jfld < gradfields.size(); ++jfld) {
        const Field& f = gradfields[jfld];
        ATLAS_ASSERT(f.functionspace() == 0 || functionspace::StructuredColumns(f.functionspace()));
    }

    const int trans_sp_nfld = compute_nfld(spfields);

    if (nfld != trans_sp_nfld) {
        throw_Exception("__invtrans_grad_adj: different number of gridpoint fields than spectral fields", Here());
    }
    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackStructuredColumns pack(rgpview);
        for (idx_t jfld = 0; jfld < gradfields.size(); ++jfld) {
            pack(gp, gradfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::InvTransAdj_t transform = ::new_invtrans_adj(trans_.get());
        transform.nscalar                = nfld;
        transform.rgp                    = rgp.data();
        transform.rspscalar              = rsp.data();
        transform.lscalarders            = true;

        TRANS_CHECK(::trans_invtrans_adj(&transform));
    }

    // Unpack the spectral fields
    {
        UnpackSpectral unpack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            unpack(spfields[jfld]);
        }
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

//---------------------------------------------------------------------------------------------

void TransIFS::__invtrans_grad_adj(const Spectral& sp, FieldSet& spfields, const functionspace::NodeColumns& gp,
                                   const FieldSet& gradfields, const eckit::Configuration& config) const {
#if ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ)
    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const idx_t nfld = compute_nfld(gradfields);
    for (idx_t jfld = 0; jfld < gradfields.size(); ++jfld) {
        const Field& f = gradfields[jfld];
        ATLAS_ASSERT(f.functionspace() == 0 || functionspace::StructuredColumns(f.functionspace()));
    }

    const int trans_sp_nfld = compute_nfld(spfields);

    if (nfld != trans_sp_nfld) {
        throw_Exception("dirtrans: different number of gridpoint fields than spectral fields", Here());
    }
    // Arrays Trans expects
    std::vector<double> rgp(nfld * ngptot());
    std::vector<double> rsp(nspec2() * nfld);
    auto rgpview = LocalView<double, 2>(rgp.data(), make_shape(nfld, ngptot()));
    auto rspview = LocalView<double, 2>(rsp.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackNodeColumns pack(rgpview, gp);
        for (idx_t jfld = 0; jfld < gradfields.size(); ++jfld) {
            pack(gradfields[jfld]);
        }
    }

    // Do transform
    {
        struct ::InvTransAdj_t transform = ::new_invtrans_adj(trans_.get());
        transform.nscalar                = nfld;
        transform.rgp                    = rgp.data();
        transform.rspscalar              = rsp.data();
        transform.lscalarders            = true;

        TRANS_CHECK(::trans_invtrans_adj(&transform));
    }

    // Unpack the spectral fields
    {
        UnpackSpectral unpack(rspview);
        for (idx_t jfld = 0; jfld < spfields.size(); ++jfld) {
            unpack(spfields[jfld]);
        }
    }
#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

// --------------------------------------------------------------------------------------------

void TransIFS::__invtrans_vordiv2wind_adj(const Spectral& sp, Field& spvor, Field& spdiv,
                                          const functionspace::StructuredColumns& gp, const Field& gpwind,
                                          const eckit::Configuration&) const {
#if ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ)

    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const size_t nfld = compute_nfld(spvor);
    if (spdiv.shape(0) != spvor.shape(0)) {
        throw_Exception("invtrans_vordiv2wind_adj: vorticity not compatible with divergence.", Here());
    }
    if (spdiv.shape(1) != spvor.shape(1)) {
        throw_Exception("invtrans_vordiv2wind_adj: vorticity not compatible with divergence.", Here());
    }
    const size_t nwindfld = compute_nfld(gpwind);
    if (nwindfld != 2 * nfld && nwindfld != 3 * nfld) {
        throw_Exception("dirtrans: wind field is not compatible with vorticity, divergence.", Here());
    }

    if (spdiv.shape(0) != nspec2()) {
        std::stringstream msg;
        msg << "invtrans_vordiv2wind_adj: Spectral vorticity and divergence have wrong dimension: "
               "nspec2 "
            << spdiv.shape(0) << " should be " << nspec2();
        throw_Exception(msg.str(), Here());
    }

    if (spvor.size() == 0) {
        throw_Exception("invtrans_vordiv2wind_adj: spectral vorticity field is empty.");
    }
    if (spdiv.size() == 0) {
        throw_Exception("invtrans_vordiv2wind_adj: spectral divergence field is empty.");
    }

    // Arrays Trans expects
    std::vector<double> rgp(2 * nfld * ngptot());
    std::vector<double> rspvor(nspec2() * nfld);
    std::vector<double> rspdiv(nspec2() * nfld);
    auto rgpview    = LocalView<double, 2>(rgp.data(), make_shape(2 * nfld, ngptot()));
    auto rspvorview = LocalView<double, 2>(rspvor.data(), make_shape(nspec2(), nfld));
    auto rspdivview = LocalView<double, 2>(rspdiv.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackStructuredColumns pack(rgpview);
        int wind_components = 2;
        pack(gp, gpwind, wind_components);
    }

    // Do transform
    {
        struct ::InvTransAdj_t transform = ::new_invtrans_adj(trans_.get());
        transform.nvordiv                = int(nfld);
        transform.rgp                    = rgp.data();
        transform.rspvor                 = rspvor.data();
        transform.rspdiv                 = rspdiv.data();

        ATLAS_ASSERT(transform.rspvor);
        ATLAS_ASSERT(transform.rspdiv);
        TRANS_CHECK(::trans_invtrans_adj(&transform));
    }

    // Pack spectral fields
    UnpackSpectral unpack_vor(rspvorview);
    unpack_vor(spvor);
    UnpackSpectral unpack_div(rspdivview);
    unpack_div(spdiv);

#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

void TransIFS::__invtrans_vordiv2wind_adj(const Spectral& sp, Field& spvor, Field& spdiv,
                                          const functionspace::NodeColumns& gp, const Field& gpwind,
                                          const eckit::Configuration&) const {
#if ATLAS_HAVE_ECTRANS || defined(TRANS_HAVE_INVTRANS_ADJ)

    assertCompatibleDistributions(gp, sp);

    // Count total number of fields and do sanity checks
    const size_t nfld = compute_nfld(spvor);
    if (spdiv.shape(0) != spvor.shape(0)) {
        throw_Exception("invtrans_vordiv2wind_adj: vorticity not compatible with divergence.", Here());
    }
    if (spdiv.shape(1) != spvor.shape(1)) {
        throw_Exception("invtrans_vordiv2wind_adj: vorticity not compatible with divergence.", Here());
    }
    const size_t nwindfld = compute_nfld(gpwind);
    if (nwindfld != 2 * nfld && nwindfld != 3 * nfld) {
        throw_Exception("dirtrans: wind field is not compatible with vorticity, divergence.", Here());
    }

    if (spdiv.shape(0) != nspec2()) {
        std::stringstream msg;
        msg << "invtrans_vordiv2wind_adj: Spectral vorticity and divergence have wrong dimension: "
               "nspec2 "
            << spdiv.shape(0) << " should be " << nspec2();
        throw_Exception(msg.str(), Here());
    }

    if (spvor.size() == 0) {
        throw_Exception("invtrans_vordiv2wind_adj: spectral vorticity field is empty.");
    }
    if (spdiv.size() == 0) {
        throw_Exception("invtrans_vordiv2wind_adj: spectral divergence field is empty.");
    }

    // Arrays Trans expects
    std::vector<double> rgp(2 * nfld * ngptot());
    std::vector<double> rspvor(nspec2() * nfld);
    std::vector<double> rspdiv(nspec2() * nfld);
    auto rgpview    = LocalView<double, 2>(rgp.data(), make_shape(2 * nfld, ngptot()));
    auto rspvorview = LocalView<double, 2>(rspvor.data(), make_shape(nspec2(), nfld));
    auto rspdivview = LocalView<double, 2>(rspdiv.data(), make_shape(nspec2(), nfld));

    // Pack gridpoints
    {
        PackNodeColumns pack(rgpview, gp);
        int wind_components = 2;
        pack(gpwind, wind_components);
    }

    // Do transform
    {
        struct ::InvTransAdj_t transform = ::new_invtrans_adj(trans_.get());
        transform.nvordiv                = int(nfld);
        transform.rgp                    = rgp.data();
        transform.rspvor                 = rspvor.data();
        transform.rspdiv                 = rspdiv.data();

        ATLAS_ASSERT(transform.rspvor);
        ATLAS_ASSERT(transform.rspdiv);
        TRANS_CHECK(::trans_invtrans_adj(&transform));
    }

    // Pack spectral fields
    UnpackSpectral unpack_vor(rspvorview);
    unpack_vor(spvor);
    UnpackSpectral unpack_div(rspdivview);
    unpack_div(spdiv);

#else
    ATLAS_NOTIMPLEMENTED;
#endif
}

///////////////////////////////////////////////////////////////////////////////

void TransIFS::distspec(const int nb_fields, const int origin[], const double global_spectra[],
                        double spectra[]) const {
    struct ::DistSpec_t args = new_distspec(trans_.get());
    args.nfld                = nb_fields;
    args.rspecg              = global_spectra;
    args.nfrom               = origin;
    args.rspec               = spectra;
    TRANS_CHECK(::trans_distspec(&args));
}

/////////////////////////////////////////////////////////////////////////////

void TransIFS::gathspec(const int nb_fields, const int destination[], const double spectra[],
                        double global_spectra[]) const {
    struct ::GathSpec_t args = new_gathspec(trans_.get());
    args.nfld                = nb_fields;
    args.rspecg              = global_spectra;
    args.nto                 = destination;
    args.rspec               = spectra;
    TRANS_CHECK(::trans_gathspec(&args));
}

/////////////////////////////////////////////////////////////////////////////

void TransIFS::distgrid(const int nb_fields, const int origin[], const double global_fields[], double fields[]) const {
    struct ::DistGrid_t args = new_distgrid(trans_.get());
    args.nfld                = nb_fields;
    args.nfrom               = origin;
    args.rgpg                = global_fields;
    args.rgp                 = fields;
    TRANS_CHECK(::trans_distgrid(&args));
}

/////////////////////////////////////////////////////////////////////////////

void TransIFS::gathgrid(const int nb_fields, const int destination[], const double fields[],
                        double global_fields[]) const {
    struct ::GathGrid_t args = new_gathgrid(trans_.get());
    args.nfld                = nb_fields;
    args.nto                 = destination;
    args.rgp                 = fields;
    args.rgpg                = global_fields;
    TRANS_CHECK(::trans_gathgrid(&args));
}

///////////////////////////////////////////////////////////////////////////////

void TransIFS::specnorm(const int nb_fields, const double spectra[], double norms[], int rank) const {
    struct ::SpecNorm_t args = new_specnorm(trans_.get());
    args.nfld                = nb_fields;
    args.rspec               = spectra;
    args.rnorm               = norms;
    args.nmaster             = rank + 1;
    TRANS_CHECK(::trans_specnorm(&args));
}

///////////////////////////////////////////////////////////////////////////////

}  // namespace trans
}  // namespace atlas
