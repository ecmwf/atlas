/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/log/Bytes.h"

#include "atlas/functionspace/Spectral.h"
#include "atlas/grid.h"
#include "atlas/linalg/dense.h"
#include "atlas/option.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/trans/LegendreCacheCreator.h"
#include "atlas/trans/Trans.h"
#include "atlas/util/Config.h"

//------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::grid;
using atlas::util::Config;
using eckit::PathName;

//------------------------------------------------------------------------------

class Tool : public AtlasTool {
    int execute(const Args& args) override;
    std::string briefDescription() override {
        return "Tool to generate a python script that plots the grid-distribution "
               "of a given grid";
    }
    std::string usage() override { return name() + " --grid=name [OPTION]... OUTPUT [--help]"; }

public:
    Tool(int argc, char** argv);

private:
    std::string key;
    PathName path_in;
    PathName path_out;
};

//-----------------------------------------------------------------------------

Tool::Tool(int argc, char** argv): AtlasTool(argc, argv) {
    add_option(new SimpleOption<std::string>(
        "grid", "Grid unique identifier\n" + indent() + "     Example values: N80, F40, O24, L32"));
    add_option(new SimpleOption<long>("truncation", "truncation (default=N-1 with N the Gaussion number)"));
    add_option(new SimpleOption<long>("nscalar", "number of scalar fields (default=1)"));
    add_option(new SimpleOption<long>("nvordiv", "number of vor/div fields (default=0)"));
    add_option(
        new SimpleOption<std::string>("type", "Type of trans implementation (if not specified, all types are used)"));
    add_option(new SimpleOption<std::string>("matrix_multiply", "backend to use in local trans type"));
    add_option(new SimpleOption<std::string>("fft", "backend to use in local trans type"));
    add_option(new SimpleOption<bool>("caching", "caching"));
    add_option(new SimpleOption<long>("niter", "number of iterations"));
}

//-----------------------------------------------------------------------------

int Tool::execute(const Args& args) {
    Trace timer(Here(), displayName());
    key = "";
    args.get("grid", key);

    std::string path_in_str = "";
    if (args.get("grid", path_in_str)) {
        path_in = path_in_str;
    }

    StructuredGrid grid;
    if (key.size()) {
        try {
            grid = Grid(key);
        }
        catch (eckit::Exception&) {
            return failed();
        }
    }
    else {
        Log::error() << "No grid specified." << std::endl;
    }

    if (!grid) {
        return failed();
    }


    int truncation;
    if (not args.get("truncation", truncation)) {
        GaussianGrid gaussian{grid};
        if (gaussian) {
            truncation = gaussian.N() - 1;  // cubic relation
        }
        else {
            Log::error() << "Could not determine truncation from grid" << std::endl;
            return failed();
        }
    }

    auto domain = GlobalDomain{};
    //auto domain = RectangularDomain{{0,90.},{-45.,45.}};


    std::vector<std::string> types;
    if (args.has("type")) {
        types.emplace_back(args.getString("type"));
        if (not trans::Trans::hasBackend(types[0])) {
            Log::error() << "Atlas was not compiled with support for '" << types[0] << "' backend." << std::endl;
            return failed();
        }
        if (not domain.global() && types[0] == "ectrans") {
            Log::error() << "The 'ectrans' backend for Trans can only use global domains." << std::endl;
            return failed();
        }
    }
    else {
        types.emplace_back("local");
        if (trans::Trans::hasBackend("ectrans") && domain.global()) {
            types.emplace_back("ectrans");
        }
    }

    bool caching  = false;
    int nb_scalar = 1;
    int nb_vordiv = 0;
    int niter     = 1;
    args.get("nscalar", nb_scalar);
    args.get("nvordiv", nb_vordiv);
    args.get("niter", niter);
    args.get("caching", caching);
    int nb_all = nb_scalar + 2 * nb_vordiv;


    Log::info() << "Configuration" << std::endl;
    Log::info() << "~~~~~~~~~~~~~" << std::endl;
    Log::info() << "  Grid           : " << grid.name() << std::endl;
    Log::info() << "  OpenMP         : " << atlas_omp_get_max_threads() << std::endl;
    Log::info() << "  Trunctation    : " << truncation << std::endl;
    Log::info() << "  trans.type     : " << types << std::endl;
    Log::info() << "  scalar fields  : " << nb_scalar << std::endl;
    Log::info() << "  vor/div fields : " << nb_vordiv << std::endl;
    Log::info() << "  niter          : " << niter << std::endl;
    Log::info() << "  caching        : " << std::boolalpha << caching << std::endl;
    if (caching) {
        Log::info() << "  cache path     : " << atlas::Library::instance().cachePath() << std::endl;
    }
    Log::info() << std::endl;

    if (caching) {
        if (not eckit::PathName{atlas::Library::instance().cachePath()}.exists()) {
            Log::error() << "Cache path " << atlas::Library::instance().cachePath() << " does not exist." << std::endl;
            Log::error() << "Set ATLAS_CACHE_PATH environment variable to an existing path." << std::endl;
            return failed();
        }
    }

    std::map<std::string, std::vector<std::string>> linalg_backends{
        {"ectrans", {"lapack"}}, {"local", {"generic", "openmp", "lapack", "eigen"}}};
    if (args.has("matrix_multiply")) {
        std::string backend      = args.getString("matrix_multiply");
        linalg_backends["local"] = {backend};
        if (not linalg::dense::Backend{backend}.available()) {
            Log::error() << "Atlas does not have support for matrix_multiply backend '" << backend << "'" << std::endl;
            return failed();
        }
    }

    std::map<std::string, std::vector<std::string>> fft_backends{
        {"ectrans", {"FFTW"}}, {"local", {"FFTW", "pocketfft"}}};
    if (args.has("fft")) {
        fft_backends["local"] = std::vector<std::string>{args.getString("fft")};
    }

    functionspace::Spectral spectral{truncation};

    std::vector<double> sp_scalar(spectral.nb_spectral_coefficients_global() * nb_scalar);
    std::vector<double> sp_vorticity(spectral.nb_spectral_coefficients_global() * nb_vordiv);
    std::vector<double> sp_divergence(spectral.nb_spectral_coefficients_global() * nb_vordiv);
    std::vector<double> gp(Grid{grid, domain}.size() * nb_all);

    for (size_t i = 0; i < nb_scalar; ++i) {
        sp_scalar[i] = 1.;
    }
    for (size_t i = 0; i < nb_vordiv; ++i) {
        sp_vorticity[i]  = 1.;
        sp_divergence[i] = 1.;
    }

    for (auto& type : types) {
        ATLAS_TRACE(type);

        trans::Cache cache;

        if (args.getBool("caching", false)) {
            trans::LegendreCacheCreator cache_creator(grid, truncation, option::type(type));
            if (cache_creator.supported()) {
                auto cachefile =
                    eckit::PathName(atlas::Library::instance().cachePath() + "/leg_" + cache_creator.uid() + ".bin");
                if (not cachefile.exists()) {
                    ATLAS_TRACE("Creating Cache");
                    Log::debug() << "Creating cache: " << cachefile
                                 << " estimated size: " << eckit::Bytes(cache_creator.estimate()) << std::endl;
                    cache_creator.create(cachefile);
                }
                Log::debug() << "Reading cache " << cachefile << " size: " << eckit::Bytes(cachefile.size())
                             << std::endl;
                ATLAS_TRACE("Reading Cache");
                cache = trans::LegendreCache(cachefile);
            }
        }

        for( auto& fft: fft_backends.at(type) ) {
            ATLAS_TRACE("fft="+fft);

            auto cstart = std::chrono::system_clock::now();
            trans::Trans trans(cache, grid, domain, truncation, option::type(type) | option::fft(fft));
            auto cend   = std::chrono::system_clock::now();
            std::chrono::duration<double> celapsed_seconds = cend - cstart;
            Log::info() << "type=" << std::setw(8) << std::left << type;
            Log::info() << "                    fft=" << std::setw(10) << std::left << fft;
            Log::info() << "  constructor:   " << celapsed_seconds.count() << " s" << std::endl;

            for (auto backend : linalg_backends.at(type)) {
                linalg::dense::current_backend(backend);
                if (linalg::dense::current_backend().available()) {
                    backend      = linalg::dense::current_backend().type();
                    auto min     = std::numeric_limits<double>::max();
                    auto max     = 0.;
                    auto zeropad = [](int n) {
                        std::stringstream s;
                        s << std::setw(3) << std::setfill('0') << n;
                        return s.str();
                    };
                    auto print = [&](const std::string& idx, double seconds) {
                        Log::info() << "type=" << std::setw(8) << std::left << type;
                        Log::info() << "  backend=" << std::setw(8) << std::left << backend;
                        Log::info() << "  fft=" << std::setw(10) << std::left << fft;
                        Log::info() << "  invtrans[" << idx << "]: " << seconds << " s" << std::endl;
                    };
                    for (size_t n = 0; n < niter; ++n) {
                        ATLAS_TRACE("invtrans [backend=" + backend + ", fft="+fft+"]");
                        auto start = std::chrono::system_clock::now();
                        trans.invtrans(nb_scalar, sp_scalar.data(), nb_vordiv, sp_vorticity.data(), sp_divergence.data(),
                                    gp.data());
                        auto end                                      = std::chrono::system_clock::now();  //
                        std::chrono::duration<double> elapsed_seconds = end - start;
                        print(zeropad(n), elapsed_seconds.count());
                        min = std::min(min, elapsed_seconds.count());
                        max = std::max(max, elapsed_seconds.count());
                    }
                    print("min", min);
                    print("max", max);
                }
            }
        }
    }

    timer.stop();
    Log::info() << Trace::report() << std::endl;
    return success();
}

//------------------------------------------------------------------------------

int main(int argc, char** argv) {
    Tool tool(argc, argv);
    return tool.start();
}
