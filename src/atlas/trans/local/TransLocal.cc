/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/trans/local/TransLocal.h"

#include <cmath>
#include <cstdlib>
#include <fstream>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/eckit.h"
#include "eckit/io/DataHandle.h"
#include "eckit/linalg/LinearAlgebra.h"
#include "eckit/log/Bytes.h"
#include "eckit/log/JSON.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library.h"
#include "atlas/linalg/dense.h"
#include "atlas/linalg/fft.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/VorDivToUV.h"
#include "atlas/trans/detail/TransFactory.h"
#include "atlas/trans/local/LegendrePolynomials.h"
#include "atlas/util/Constants.h"
#include "pluto/pluto.h"

#include "atlas/library/defines.h"
#if ATLAS_HAVE_FFTW
#include "fftw3.h"
#endif

// move latitudes at the poles to the following latitude:
// (otherwise we would divide by zero when computing u,v from U,V)
static constexpr double latPole = 89.9999999;
// (latPole=89.9999999 seems to produce the best accuracy. Moving it further away
//  or closer to the pole both increase the errors!)

namespace atlas {
namespace trans {

namespace {
static TransBuilderGrid<TransLocal> builder("local", "local");
}  // namespace

namespace {
class TransParameters {
public:
    TransParameters(const eckit::Configuration& config): config_(config) {}
    ~TransParameters() = default;

    /*
     * For the future
     */
    //    bool scalar_derivatives() const { return config_.getBool( "scalar_derivatives", false ); }

    //    bool wind_EW_derivatives() const { return config_.getBool( "wind_EW_derivatives", false ); }

    //    bool vorticity_divergence_fields() const { return config_.getBool( "vorticity_divergence_fields", false ); }

    //    std::string read_legendre() const { return config_.getString( "read_legendre", "" ); }
    //    bool global() const { return config_.getBool( "global", false ); }

    //    std::string read_fft() const { return config_.getString( "read_fft", "" ); }


    std::string write_legendre() const { return config_.getString("write_legendre", ""); }

    std::string write_fft() const { return config_.getString("write_fft", ""); }


    bool export_legendre() const { return config_.getBool("export_legendre", false); }

    int warning() const { return config_.getInt("warning", 1); }

    std::string fft() const {
        static const std::set<std::string> supported = {"OFF","FFTW","pocketfft"};
        std::string defaultLinalgFFTBackend = atlas::Library::instance().linalgFFTBackend();
        if (defaultLinalgFFTBackend.empty()) {
            defaultLinalgFFTBackend = "FFTW";
        }
        std::string value = config_.getString("fft", defaultLinalgFFTBackend);
        if (supported.find(value) == supported.end()) {
            ATLAS_THROW_EXCEPTION("FFT backend \"" << value << "\" is not one of the supported : OFF, FFTW, pocketfft");
        }
        if (not ATLAS_HAVE_FFTW && value == "FFTW") {
            value = "pocketfft";
        }
        return value;
    }

    std::string matrix_multiply() const { return config_.getString("matrix_multiply", ""); }


private:
    const eckit::Configuration& config_;
};

struct ReadCache {
    ReadCache(const void* cache) {
        begin = reinterpret_cast<const char*>(cache);
        pos   = 0;
    }
    template <typename T>
    T* read(size_t size) {
        const T* v = reinterpret_cast<const T*>(begin + pos);
        pos += size * sizeof(T);
        return const_cast<T*>(v);
    }

    const char* begin;
    size_t pos;
};

struct WriteCache {
    WriteCache(const eckit::PathName& file_path): dh_(file_path.fileHandle(/*overwrite = */ true)) {
        if (file_path.exists()) {
            std::stringstream err;
            err << "Cannot open cache file " << file_path << " for writing as it already exists. Remove first.";
            throw_Exception(err.str(), Here());
        }
        dh_->openForWrite(0);
        pos = 0;
    }
    ~WriteCache() { dh_->close(); }
    template <typename T>
    void write(const T* v, long size) {
        dh_->write(v, size * sizeof(T));
        pos += size * sizeof(T);
    }

    //void write( long v ) {
    //    dh_->write( &v , sizeof(long) );
    //    pos += sizeof(long);
    //}

    //void write( const Grid& grid ) {
    //    std::stringstream s;
    //    eckit::JSON json(s);
    //    json << grid.spec();
    //    std::string grid_spec( s.str() );
    //    long size = grid_spec.size();
    //    write( size );
    //    dh_->write( grid_spec.c_str(), grid_spec.size() );
    //    pos += grid_spec.size();
    //}


    std::unique_ptr<eckit::DataHandle> dh_;
    size_t pos;
};

}  // namespace

// --------------------------------------------------------------------------------------------------------------------
// Helper functions
// --------------------------------------------------------------------------------------------------------------------
namespace {  // anonymous

size_t legendre_size(const size_t truncation) {
    return (truncation + 2) * (truncation + 1) / 2;
}

//int nlats_northernHemisphere( const int nlats ) {
//    return ceil( nlats / 2. );
//    // using ceil here should make it possible to have odd number of latitudes (with the centre latitude being the equator)
//}

size_t num_n(const int truncation, const int m, const bool symmetric) {
    int len = (truncation - m + (symmetric ? 2 : 1)) / 2;
    ATLAS_ASSERT(len >= 0);
    return size_t(len);
}

size_t add_padding(size_t n) {
    return size_t(std::ceil(n / 8.)) * 8;
}

std::string detect_linalg_backend(const std::string& linalg_backend_) {
    linalg::dense::Backend linalg_backend = linalg::dense::Backend{linalg_backend_};
    if (linalg_backend.type() == linalg::dense::backend::eckit_linalg::type()) {
        std::string backend;
        linalg_backend.get("backend", backend);
        if (backend.empty() || backend == "default") {
#if ATLAS_ECKIT_HAVE_ECKIT_585
            return eckit::linalg::LinearAlgebraDense::backend().name();
#else
            return eckit::linalg::LinearAlgebra::backend().name();
#endif
        }
        return backend;
    }
    return linalg_backend.type();
};

bool using_eckit_default_backend(const std::string& linalg_backend_) {
    linalg::dense::Backend linalg_backend = linalg::dense::Backend{linalg_backend_};
    if (linalg_backend.type() == linalg::dense::backend::eckit_linalg::type()) {
        std::string backend;
        linalg_backend.get("backend", backend);
        if (backend.empty() || backend == "default") {
            return true;
        }
    }
    return false;
};

static constexpr size_t alignment = 256;

template <typename T>
static void allocate_aligned(pluto::memory_resource* mr, T*& ptr, size_t n, std::string_view label = "") {
    ptr = pluto::allocator<T>{mr}.allocate(label, n);
}

template <typename T>
static void deallocate_aligned(pluto::memory_resource* mr, T*& ptr, size_t n, std::string_view label = "") {
    pluto::allocator<T>{mr}.deallocate(label, ptr, n);
    ptr = nullptr;
}

}  // namespace

pluto::memory_resource* TransLocal::memory_resource() const {
    // return memory_resource_.get();
    return pluto::host_pool_resource();
}

void TransLocal::alloc_aligned(double*& ptr, size_t n, const char* msg, pluto::memory_resource* mr) const {
    ATLAS_ASSERT(msg);
    Log::debug() << "TransLocal: allocating '" << msg << "': " << eckit::Bytes(sizeof(double) * n) << std::endl;
    if (mr == nullptr) {
        mr = memory_resource();
    }
    allocate_aligned(mr, ptr, n, msg);
}

void TransLocal::free_aligned(double*& ptr, size_t n, const char* msg, pluto::memory_resource* mr) const {
    ATLAS_ASSERT(msg);
    Log::debug() << "TransLocal: deallocating '" << msg << "'" << std::endl;
    if (mr == nullptr) {
        mr = memory_resource();
    }
    deallocate_aligned(mr, ptr, n, msg);
}

int fourier_truncation(const int truncation,   // truncation
                       const int nx,           // number of longitudes
                       const int /*nxmax*/,    // maximum nx
                       const int ndgl,         // number of latitudes
                       const double lat,       // latitude in radian
                       const bool fullgrid) {  // regular grid
    int trc     = truncation;
    int trclin  = ndgl - 1;
    int trcquad = ndgl * 2 / 3 - 1;
    if (truncation >= trclin || fullgrid) {
        // linear
        trc = (nx - 1) / 2;
    }
    else if (truncation >= trcquad) {
        // quadratic
        double weight = 3 * (trclin - truncation) / ndgl;
        double sqcos  = std::pow(std::cos(lat), 2);

        trc = static_cast<int>((nx - 1) / (2 + weight * sqcos));
    }
    else {
        // cubic
        double sqcos = std::pow(std::cos(lat), 2);

        trc = static_cast<int>((nx - 1) / (2 + sqcos) - 1);
    }
    trc = std::min(truncation, trc);
    return trc;
}

namespace detail {

struct FFT_Data {
    mdspan<std::complex<double>,dims<2>,layout_stride> many_in() {
        return mdspan<std::complex<double>,dims<2>,layout_stride>{in_, {dims<2>{howmany_,size_in_},std::array{dist_in_,std::size_t{1}}}};
    }
    mdspan<double,dims<2>,layout_stride> many_out() {
        return mdspan<double,dims<2>,layout_stride>{out_, {dims<2>{howmany_,size_out_},std::array{dist_out_,std::size_t{1}}}};
    }

    std::complex<double>* one_in(size_t j = 0) {
        ATLAS_ASSERT(j < howmany_);
        return in_ + j * dist_in_;
    }
    double* one_out(size_t j = 0) {
        ATLAS_ASSERT(j < howmany_);
        return out_ + j * dist_out_;
    }

    void allocate(size_t size_out, size_t howmany, pluto::memory_resource* mr = nullptr) {
        howmany_  = howmany;
        size_out_ = size_out;
        size_in_  = (size_out_ / 2) + 1;
        dist_in_  = pluto::align_up(size_in_,  alignment/sizeof(std::complex<double>));
        dist_out_ = pluto::align_up(size_out_, alignment/sizeof(double));
        if (mr == nullptr) {
            mr = pluto::host::get_default_resource();
        }
        mr_ = mr;
        allocate_aligned(mr_, in_,  howmany_*dist_in_, "fft_data_.in_");
        allocate_aligned(mr_, out_, howmany_*dist_out_, "fft_data_.out");
    }
    void deallocate() {
        deallocate_aligned(mr_, in_,  howmany_ * dist_in_);
        deallocate_aligned(mr_, out_, howmany_ * dist_out_);
    }
private:
    std::complex<double>* in_{nullptr}; // input of inverse FFT, complex
    double* out_{nullptr}; // output of inverse FFT, real
    std::size_t size_in_;  // size of 1 FFT, complex in
    std::size_t size_out_; // size of 1 FFT, real out
    std::size_t dist_in_;  // padded size between consecutive FFTs, complex in
    std::size_t dist_out_; // padded size between consecutive FFTs, real out
    std::size_t howmany_;  // howmany FFTs are allocated with padding
    pluto::memory_resource* mr_;
};

}  // namespace detail


// --------------------------------------------------------------------------------------------------------------------
// Class TransLocal
// --------------------------------------------------------------------------------------------------------------------

bool TransLocal::warning(const eckit::Configuration& config) const {
    int warning = warning_;
    config.get("warning", warning);
    return (warning > 0 && grid_.size() >= warning);
}

TransLocal::TransLocal(const Cache& cache, const Grid& grid, const Domain& domain, const long truncation,
                       const eckit::Configuration& config):
    grid_(grid, domain),
    truncation_(static_cast<int>(truncation)),
    precompute_(config.getBool("precompute", true)),
    cache_(cache),
    legendre_cache_(cache.legendre().data()),
    legendre_cachesize_(cache.legendre().size()),
    fft_cache_(cache.fft().data()),
    fft_cachesize_(cache.fft().size()),
    fft_data_(new detail::FFT_Data),
    linalg_backend_(TransParameters{config}.matrix_multiply()),
    fft_backend_(TransParameters{config}.fft()),
    warning_(TransParameters{config}.warning()),
    upstream_memory_resource_(pluto::host::get_default_resource()) {
    ATLAS_TRACE("TransLocal constructor");

    if (mpi::size() > 1) {
        ATLAS_THROW_EXCEPTION("TransLocal is not implemented for more than 1 MPI task.");
    }

    double fft_threshold = 0.0;  // fraction of latitudes of the full grid down to which FFT is used.
    // This threshold needs to be adjusted depending on the dgemm and FFT performance of the machine
    // on which this code is running!
    idx_t nlats       = 0;
    idx_t nlonsMax    = 0;
    idx_t neqtr       = 0;
    useFFT_           = (fft_backend_ != "OFF");
    if (useFFT_) {
        if (fft_backend_ == "FFTW") {
            fft_.reset(new linalg::FFTW());
        }
        else if(fft_backend_ == "pocketfft") {
            fft_.reset(new linalg::pocketfft());
        }
        else {
            ATLAS_THROW_EXCEPTION("Unrecognised FFT configuration");
        }
    }

    unstruct_precomp_ = (config.has("precompute") ? precompute_ : false);
    no_symmetry_      = false;
    nlatsNH_          = 0;
    nlatsSH_          = 0;
    nlatsLeg_         = 0;
    nlatsLegDomain_   = 0;
    nlatsLegReduced_  = 0;
    bool useGlobalLeg = true;
    bool no_nest      = false;

    if (StructuredGrid(grid_) && not grid_.projection()) {
        StructuredGrid g(grid_);
        nlats    = g.ny();
        nlonsMax = g.nxmax();

        // check location of domain relative to the equator:
        for (idx_t j = 0; j < nlats; ++j) {
            // assumptions: latitudes in g.y(j) are monotone and decreasing
            // no assumption on whether we have 0, 1 or 2 latitudes at the equator
            double lat = g.y(j);
            (eckit::types::is_approximately_equal(lat, 0.) ? neqtr : (lat < 0 ? nlatsSH_ : nlatsNH_))++;
        }
        if (neqtr > 0) {
            nlatsNH_++;
            nlatsSH_++;
        }
        if (nlatsNH_ >= nlatsSH_) {
            nlatsLegDomain_ = nlatsNH_;
        }
        else {
            nlatsLegDomain_ = nlatsSH_;
        }

        gridGlobal_ = grid;
        if (not gridGlobal_.domain().global()) {
            // The grid is not a nest of a global grid
            if (RegularGrid(grid_)) {
                no_nest         = true;
                no_symmetry_    = true;
                useFFT_         = false;
                fft_backend_    = "OFF";
                nlatsNH_        = nlats;
                nlatsSH_        = 0;
                nlatsLegDomain_ = nlatsNH_;
                gridGlobal_     = grid_;
                useGlobalLeg    = false;
            }
            else {  // non-nested reduced grids are not supported
                std::ostringstream log;
                log << "Transforms to non-regular regional grids is not supported, unless it defined as a cropping of "
                       "a global grid."
                    << std::endl;
                log << "    Input grid: " << grid.spec() << std::endl;
                log << "    Applied domain: " << domain << std::endl;
                log << "    Regional grid after Domain applied (the output grid): " << grid_.spec() << std::endl;
                throw_NotImplemented(log.str(), Here());
            }
        }

        StructuredGrid gs_global(gridGlobal_);
        ATLAS_ASSERT(gs_global);  // assert structured grid
        StructuredGrid gsLeg = (useGlobalLeg ? gs_global : g);
        nlonsMaxGlobal_      = gs_global.nxmax();
        jlonMin_.resize(1);
        jlonMin_[0]  = 0;
        jlatMin_     = 0;
        nlatsGlobal_ = gs_global.ny();
        if (grid_.domain().global()) {
            Log::debug() << "Global grid with " << nlats << " latitudes." << std::endl;
        }
        else {
            Log::debug() << "Grid has " << nlats << " latitudes. Global grid has " << nlatsGlobal_ << std::endl;
        }
        if (useGlobalLeg) {
            nlatsLeg_ = (nlatsGlobal_ + 1) / 2;
        }
        else {
            nlatsLeg_        = nlatsLegDomain_;
            nlatsLegReduced_ = nlatsLeg_;
        }
        for (int jlat = 0; jlat < nlatsGlobal_; jlat++) {
            if (gs_global.y(jlat) > g.y(0)) {
                //Log::info() << gs_global.y( jlat ) << ">" << g.y( 0 ) << " ";
                jlatMin_++;
            };
        }
        //Log::info() << std::endl;
        int jlatMinLeg_ = jlatMin_;
        if (nlatsNH_ < nlatsSH_) {
            jlatMinLeg_ += nlatsNH_ - nlatsSH_;
        };
        if (jlatMin_ >= (nlatsGlobal_ + 1) / 2) {
            jlatMinLeg_ -= 2 * (jlatMin_ - (nlatsGlobal_ + 1) / 2);
            if (nlatsGlobal_ % 2 == 1) {
                jlatMinLeg_--;
            }
        };
        if (useGlobalLeg) {
            nlatsLegReduced_ = jlatMinLeg_ + nlatsLegDomain_;
        }

        // reduce truncation towards the pole for reduced meshes:
        nlat0_.resize(truncation_ + 1);
        if (no_nest) {
            for (int j = 0; j <= truncation_; j++) {
                nlat0_[j] = 0;
            }
        }
        else {
            int nmen0 = -1;
            for (int jlat = 0; jlat < nlatsGlobal_ / 2; jlat++) {
                double lat = gs_global.y(jlat) * util::Constants::degreesToRadians();
                int nmen   = fourier_truncation(truncation_, gs_global.nx(jlat), gs_global.nxmax(), nlatsGlobal_, lat,
                                                RegularGrid(gs_global));
                nmen       = std::max(nmen0, nmen);
                int ndgluj = nlatsLeg_ - std::min(nlatsLeg_, nlatsLeg_ + jlatMinLeg_ - jlat);
                if (useGlobalLeg) {
                    ndgluj = std::max(jlatMinLeg_, jlat);
                }
                for (int j = nmen0 + 1; j <= nmen; j++) {
                    nlat0_[j] = ndgluj;
                }
                nmen0 = nmen;
            }
            for (int j = nmen0 + 1; j <= truncation_; j++) {
                nlat0_[j] = nlatsLeg_;
            }
        }
        /*Log::info() << "nlats=" << g.ny() << " nlatsGlobal=" << gs_global.ny() << " jlatMin=" << jlatMin_
                    << " jlatMinLeg=" << jlatMinLeg_ << " nlatsGlobal/2-nlatsLeg=" << nlatsGlobal_ / 2 - nlatsLeg_
                    << " nlatsLeg_=" << nlatsLeg_ << " nlatsLegDomain_=" << nlatsLegDomain_ << std::endl;*/

        // compute longitudinal location of domain within global grid for using FFT:
        auto wrapAngle = [&](double angle) {
            double result = std::fmod(angle, 360.);
            if (result < 0.) {
                result += 360.;
            }
            return result;
        };
        if (useFFT_) {
            double lonmin = wrapAngle(g.x(0, 0));
            if (nlonsMax < fft_threshold * nlonsMaxGlobal_) {
                useFFT_ = false;
                fft_backend_ = "OFF";
            }
            else {
                // need to use FFT with cropped grid
                if (RegularGrid(gridGlobal_)) {
                    for (idx_t jlon = 0; jlon < nlonsMaxGlobal_; ++jlon) {
                        if (gs_global.x(jlon, 0) < lonmin) {
                            jlonMin_[0]++;
                        }
                    }
                }
                else {
                    nlonsGlobal_.resize(nlats);
                    jlonMin_.resize(nlats);
                    atlas_omp_parallel_for (idx_t jlat = 0; jlat < nlats; jlat++) {
                        double lonmin      = wrapAngle(g.x(0, jlat));
                        nlonsGlobal_[jlat] = gs_global.nx(jlat + jlatMin_);
                        jlonMin_[jlat]     = 0;
                        for (idx_t jlon = 0; jlon < nlonsGlobal_[jlat]; ++jlon) {
                            if (gs_global.x(jlon, jlat + jlatMin_) < lonmin) {
                                jlonMin_[jlat]++;
                            }
                        }
                    }
                }
            }
        }
        //Log::info() << "nlats=" << g.ny() << " nlatsGlobal=" << gs_global.ny() << std::endl;
        std::vector<double> lats(nlatsLeg_);
        std::vector<double> lons(nlonsMax);
        if (nlatsNH_ >= nlatsSH_ || useGlobalLeg) {
            atlas_omp_parallel_for (idx_t j = 0; j < nlatsLeg_; ++j) {
                double lat = gsLeg.y(j);
                if (lat > latPole) {
                    lat = latPole;
                }
                if (lat < -latPole) {
                    lat = -latPole;
                }
                lats[j] = lat * util::Constants::degreesToRadians();
            }
        }
        else {
            atlas_omp_parallel_for (idx_t idx = 0; idx < nlatsLeg_; ++idx) {
                idx_t j = nlats-1 - idx;
                double lat = gsLeg.y(j);
                if (lat > latPole) {
                    lat = latPole;
                }
                if (lat < -latPole) {
                    lat = -latPole;
                }
                lats[idx] = -lat * util::Constants::degreesToRadians();
            }
        }
        atlas_omp_parallel_for (idx_t j = 0; j < nlonsMax; ++j) {
            lons[j] = g.x(j, 0) * util::Constants::degreesToRadians();
        }
        /*Log::info() << "lats: ";
        for ( int j = 0; j < nlatsLeg_; j++ ) {
            Log::info() << lats[j] << " ";
        }
        Log::info() << std::endl;*/


        Log::debug() << "TransLocal set up with:\n"
                     << " - grid: " << grid.name() << '\n'
                     << " - truncation: " << truncation << '\n';
        if (not domain.global()) {
            Log::debug() << " - domain: " << domain << '\n';
        }
        if (GlobalDomain(domain)) {
            if (GlobalDomain(domain).west() != 0.) {
                Log::debug() << " - global domain with modified west: " << GlobalDomain(domain).west() << '\n';
            }
        }
        Log::debug() << " - fft: " << fft_backend_ << '\n';
        Log::debug() << " - linalg_backend: ";
        if (using_eckit_default_backend(linalg_backend_)) {
            Log::debug() << "eckit_linalg default (currently \"" << detect_linalg_backend(linalg_backend_)
                         << "\" but could be changed after setup, check invtrans debug output)" << '\n';
        }
        else {
            Log::debug() << detect_linalg_backend(linalg_backend_) << '\n';
        }
        Log::debug() << " - legendre_cache: " << std::boolalpha << bool(legendre_cache_) << std::endl;


        // precomputations for Legendre polynomials:
        {
            const auto nlatsLeg = size_t(nlatsLeg_);
            size_t size_sym     = 0;
            size_t size_asym    = 0;
            legendre_sym_begin_.resize(truncation_ + 3);
            legendre_asym_begin_.resize(truncation_ + 3);
            legendre_sym_begin_[0]  = 0;
            legendre_asym_begin_[0] = 0;
            for (idx_t jm = 0; jm <= truncation_ + 1; jm++) {
                size_sym  += add_padding(num_n(truncation_ + 1, jm, /*symmetric*/     true ) * nlatsLeg);
                size_asym += add_padding(num_n(truncation_ + 1, jm, /*antisymmetric*/ false) * nlatsLeg);
                legendre_sym_begin_[jm + 1]  = size_sym;
                legendre_asym_begin_[jm + 1] = size_asym;
            }

            if (legendre_cache_) {
                ReadCache legendre(legendre_cache_);
                legendre_sym_  = legendre.read<double>(size_sym);
                legendre_asym_ = legendre.read<double>(size_asym);
                ATLAS_ASSERT(legendre.pos == legendre_cachesize_);
                // TODO: check this is all aligned...
            }
            else {
                if (TransParameters(config).export_legendre()) {
                    ATLAS_ASSERT(not cache_.legendre());

                    size_t bytes = sizeof(double) * (size_sym + size_asym);
                    Log::debug() << "TransLocal: allocating LegendreCache: " << eckit::Bytes(bytes) << std::endl;
                    export_legendre_ = LegendreCache(bytes);

                    legendre_cachesize_ = export_legendre_.legendre().size();
                    legendre_cache_     = export_legendre_.legendre().data();
                    ReadCache legendre(legendre_cache_);
                    legendre_sym_  = legendre.read<double>(size_sym);
                    legendre_asym_ = legendre.read<double>(size_asym);
                }
                else {
                    legendre_sym_size_ = size_sym;
                    legendre_asym_size_ = size_asym;
                    alloc_aligned(legendre_sym_,  size_sym,  "Legendre coeffs symmetric",  upstream_memory_resource_);
                    alloc_aligned(legendre_asym_, size_asym, "Legendre coeffs asymmetric", upstream_memory_resource_);
                }

                ATLAS_TRACE_SCOPE("Legendre precomputations (structured)") {
                    compute_legendre_polynomials(truncation_ + 1, nlatsLeg_, lats.data(), legendre_sym_, legendre_asym_,
                                                 legendre_sym_begin_.data(), legendre_asym_begin_.data());
                }
                std::string file_path = TransParameters(config).write_legendre();
                if (file_path.size()) {
                    ATLAS_TRACE("Write LegendreCache to file");
                    Log::debug() << "Writing Legendre cache file ..." << std::endl;
                    Log::debug() << "    path: " << file_path << std::endl;
                    WriteCache legendre(file_path);
                    legendre.write(legendre_sym_, size_sym);
                    legendre.write(legendre_asym_, size_asym);
                    Log::debug() << "    size: " << eckit::Bytes(legendre.pos) << std::endl;
                }
            }
        }

        // precomputations for Fourier transformations:
        if (useFFT_) {
            ATLAS_TRACE("Fourier precomputations ("+fft_->type()+")");
            fft_data_->allocate(nlonsMaxGlobal_, (gridGlobal_) ? nlats : atlas_omp_get_max_threads(), upstream_memory_resource_);

#if ATLAS_HAVE_FFTW
            if (fft_cache_) {
                Log::debug() << "Import FFTW wisdom from cache" << std::endl;
                fftw_import_wisdom_from_string(static_cast<const char*>(fft_cache_));
            }
                //                std::string wisdomString( "" );
                //                std::ifstream read( "wisdom.bin" );
                //                if ( read.is_open() ) {
                //                    std::getline( read, wisdomString );
                //                    while ( read ) {
                //                        std::string line;
                //                        std::getline( read, line );
                //                        wisdomString += line;
                //                    }
                //                }
                //                read.close();
                //                if ( wisdomString.length() > 0 ) { fftw_import_wisdom_from_string( &wisdomString[0u] ); }
#endif

            if (RegularGrid(gridGlobal_)) {
                fft_->prepare_inverse_c2r_many(fft_data_->many_in(), fft_data_->many_out());
            }
            else {
                for (int j = 0; j < nlatsLegDomain_; j++) {
                    int nlonsGlobalj = gs_global.nx(jlatMinLeg_ + j);
                    // Creating plans is not thread-safe
                    fft_->prepare_inverse_c2r(nlonsGlobalj, fft_data_->one_in(), fft_data_->one_out());
                }
            }

#if ATLAS_HAVE_FFTW
            std::string file_path = TransParameters(config).write_fft();
            if (file_path.size()) {
                Log::debug() << "Write FFTW wisdom to file " << file_path << std::endl;
                //bool success = fftw_export_wisdom_to_filename( "wisdom.bin" );
                //ASSERT( success );
                //std::ofstream write( file_path );
                //write << FFTW_Wisdom();

                FILE* file_fftw = fopen(file_path.c_str(), "wb");
                fftw_export_wisdom_to_file(file_fftw);
                fclose(file_fftw);
            }
            //                std::string newWisdom( fftw_export_wisdom_to_string() );
            //                if ( 1.1 * wisdomString.length() < newWisdom.length() ) {
            //                    std::ofstream write( "wisdom.bin" );
            //                    write << newWisdom;
            //                    write.close();
            //                }
#else
            // std::string file_path = TransParameters(config).write_fft();
            // if (file_path.size()) {
            //     std::ofstream write(file_path);
            //     write << "No cache available, as FFTW is not enabled" << std::endl;
            //     write.close();
            // }
#endif
        }
        else { // not using FFT
            Log::warning()
                << "WARNING: Spectral transform results may contain aliasing errors. This will be addressed soon."
                << std::endl;
            fourier_size_ =  2 * (truncation_ + 1) * nlonsMax;
            alloc_aligned(fourier_, fourier_size_, "Fourier coeffs", upstream_memory_resource_);
            {
                ATLAS_TRACE("Fourier precomputations (NoFFT)");
                int idx = 0;
                for (int jm = 0; jm < truncation_ + 1; jm++) {
                    double factor = 1.;
                    if (jm > 0) {
                        factor = 2.;
                    }
                    for (int jlon = 0; jlon < nlonsMax; jlon++) {
                        fourier_[idx++] = +std::cos(jm * lons[jlon]) * factor;  // real part
                    }
                    for (int jlon = 0; jlon < nlonsMax; jlon++) {
                        fourier_[idx++] = -std::sin(jm * lons[jlon]) * factor;  // imaginary part
                    }
                }
            }
        }
    }
    else {
        // unstructured grid
        if (unstruct_precomp_) {
            ATLAS_TRACE("Legendre precomputations (unstructured)");

            if (warning()) {
                Log::warning()
                    << "WARNING: Precomputations for spectral transforms could take a long time as there's no structure"
                       " to take advantage of!!!"
                    << std::endl
                    << "The precomputed values consume at least "
                    << eckit::Bytes(sizeof(double) * legendre_size(truncation_) * grid_.size()) << " ("
                    << eckit::Bytes(sizeof(double) * legendre_size(truncation_)) << " for each of " << grid_.size()
                    << " grid points )" << std::endl
                    << "Furthermore, results may contain aliasing errors." << std::endl;
            }

            std::vector<double, pluto::allocator<double>> lats(grid_.size(), memory_resource());
            legendre_size_ = legendre_size(truncation_) * grid_.size();
            alloc_aligned(legendre_, legendre_size(truncation_) * grid_.size(), "Legendre coeffs.", upstream_memory_resource_);
            int j(0);
            for (PointLonLat p : grid_.lonlat()) {
                lats[j++] = p.lat() * util::Constants::degreesToRadians();
            }
            compute_legendre_polynomials_all(truncation_, grid_.size(), lats.data(), legendre_);
        }
        if (TransParameters(config).write_legendre().size()) {
            throw_NotImplemented(
                "Caching for unstructured grids or structured grids with projections not yet implemented", Here());
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

TransLocal::TransLocal(const Grid& grid, const long truncation, const eckit::Configuration& config):
    TransLocal(Cache(), grid, grid.domain(), truncation, config) {}

TransLocal::TransLocal(const Grid& grid, const Domain& domain, const long truncation,
                       const eckit::Configuration& config):
    TransLocal(Cache(), grid, domain, truncation, config) {}

TransLocal::TransLocal(const Cache& cache, const Grid& grid, const long truncation, const eckit::Configuration& config):
    TransLocal(cache, grid, grid.domain(), truncation, config) {}

// --------------------------------------------------------------------------------------------------------------------

TransLocal::~TransLocal() {
    if (StructuredGrid(grid_) && not grid_.projection()) {
        if (not legendre_cache_) {
            free_aligned(legendre_sym_,  legendre_sym_size_,  "symmetric", upstream_memory_resource_);
            free_aligned(legendre_asym_, legendre_asym_size_, "asymmetric", upstream_memory_resource_);
        }
        if (useFFT_) {
            fft_data_->deallocate();
        }
        else {
            free_aligned(fourier_, fourier_size_, "Fourier coeffs", upstream_memory_resource_);
        }
    }
    else {
        if (unstruct_precomp_) {
            free_aligned(legendre_, legendre_size_, "Legendre coeffs", upstream_memory_resource_);
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

const functionspace::Spectral& TransLocal::spectral() const {
    if (not spectral_) {
        spectral_ = functionspace::Spectral( truncation() );
    }
    return spectral_;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans(const Field& spfield, Field& gpfield, const eckit::Configuration& config) const {
    // VERY PRELIMINARY IMPLEMENTATION WITHOUT ANY GUARANTEES
    int nb_scalar_fields = 1;
    ATLAS_ASSERT(spfield.rank() == 1, "Only rank-1 fields supported at the moment");
    ATLAS_ASSERT(gpfield.rank() == 1, "Only rank-1 fields supported at the moment");
    const auto scalar_spectra = array::make_view<double, 1>(spfield);
    auto gp_fields            = array::make_view<double, 1>(gpfield);

    if (gp_fields.shape(0) < grid().size()) {
        // Hopefully the halo (if present) is appended
        ATLAS_DEBUG_VAR(gp_fields.shape(0));
        ATLAS_DEBUG_VAR(grid().size());
        ATLAS_ASSERT(gp_fields.shape(0) < grid().size());
    }

    invtrans(nb_scalar_fields, scalar_spectra.data(), gp_fields.data(), config);
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans(const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& config) const {
    // VERY PRELIMINARY IMPLEMENTATION WITHOUT ANY GUARANTEES
    ATLAS_ASSERT(spfields.size() == gpfields.size());
    for (idx_t f = 0; f < spfields.size(); ++f) {
        invtrans(spfields[f], gpfields[f], config);
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_grad(const Field& /*spfield*/, Field& /*gradfield*/, const eckit::Configuration&) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_grad(const FieldSet& /*spfields*/, FieldSet& /*gradfields*/,
                               const eckit::Configuration&) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void gp_transpose(const int nb_size, const int nb_fields, const double gp_tmp[], double gp_fields[]) {
    for (int jgp = 0; jgp < nb_size; jgp++) {
        for (int jfld = 0; jfld < nb_fields; jfld++) {
            gp_fields[jfld * nb_size + jgp] = gp_tmp[jgp * nb_fields + jfld];
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_vordiv2wind(const Field& spvor, const Field& spdiv, Field& gpwind,
                                      const eckit::Configuration& config) const {
    // VERY PRELIMINARY IMPLEMENTATION WITHOUT ANY GUARANTEES
    ATLAS_ASSERT(spvor.rank() == 1, "Only rank-1 fields supported at the moment");
    ATLAS_ASSERT(spdiv.rank() == 1, "Only rank-1 fields supported at the moment");
    int nb_vordiv_fields          = 1;
    const auto vorticity_spectra  = array::make_view<double, 1>(spvor);
    const auto divergence_spectra = array::make_view<double, 1>(spdiv);
    auto gp_fields                = array::make_view<double, 2>(gpwind);

    size_t spectral_data_size = 2 * legendre_size(truncation_) * nb_vordiv_fields;
    ATLAS_ASSERT(vorticity_spectra.size() == spectral_data_size);
    ATLAS_ASSERT(divergence_spectra.size() == spectral_data_size);

    if (gp_fields.shape(1) == grid().size() && gp_fields.shape(0) == 2) {
        invtrans(nb_vordiv_fields, vorticity_spectra.data(), divergence_spectra.data(), gp_fields.data(), config);
    }
    else if (gp_fields.shape(0) == grid().size() && gp_fields.shape(1) == 2) {
        array::ArrayT<double> gpwind_t(gp_fields.shape(1), gp_fields.shape(0));
        auto gp_fields_t = array::make_view<double, 2>(gpwind_t);
        invtrans(nb_vordiv_fields, vorticity_spectra.data(), divergence_spectra.data(), gp_fields_t.data(), config);
        gp_transpose(grid().size(), 2, gp_fields_t.data(), gp_fields.data());
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void TransLocal::invtrans_adj(const Field& gpfield, Field& spfield, const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_adj(const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_grad_adj(const Field& /*gradfield*/, Field& /*spfield*/, const eckit::Configuration&) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_grad_adj(const FieldSet& /*gradfields*/, FieldSet& /*spfields*/,
                                   const eckit::Configuration&) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_vordiv2wind_adj(const Field& gpwind, Field& spvor, Field& spdiv,
                                          const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans(const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                          const eckit::Configuration& config) const {
    invtrans_uv(truncation_, nb_scalar_fields, 0, scalar_spectra, gp_fields, config);
}


// --------------------------------------------------------------------------------------------------------------------

// This is a C++20 std function!
template <std::size_t N, typename T>
#if defined(__clang__) || defined(__GNUC__)
  __attribute__((always_inline))
#endif
constexpr T* assume_aligned(T* ptr) {
#if defined(__clang__) || (defined(__GNUC__) && !defined(__ICC))
    return reinterpret_cast<T*>(__builtin_assume_aligned(ptr, N));
#else
    return ptr;
#endif
}

void TransLocal::invtrans_legendre(const int truncation, const int nlats, const int nb_fields,
                                   const int /*nb_vordiv_fields*/, const double scalar_spectra[], double scl_fourier[],
                                   const eckit::Configuration&) const {
    // Legendre transform:
    Log::debug() << "TransLocal::invtrans_legendre: Legendre GEMM with \"" << detect_linalg_backend(linalg_backend_)
                    << "\" using " << nlatsLegReduced_ - nlat0_[0] << " latitudes out of " << nlatsGlobal_ / 2
                    << std::endl;
    linalg::dense::Backend linalg_backend{linalg_backend_};

    using fourier_out_mdspan_t = mdspan<double,extents<int,dynamic_extent,dynamic_extent,dynamic_extent,2>>;
    fourier_out_mdspan_t fourier_out{scl_fourier, nb_fields, nlats, truncation_ + 1};

    // Preallocate some reusable buffers for within the jm loop
    size_t scl_fourier_size = 0;
    size_t scalar_sym_size  = 0;
    size_t scalar_asym_size = 0;
    {
        ATLAS_TRACE("find sizes");
        for (int jm = 0; jm <= truncation_; jm++) {
            size_t size_sym  = num_n(truncation_ + 1, jm, true);
            size_t size_asym  = num_n(truncation_ + 1, jm, false);
            const int n_imag = (jm ? 2 : 1);
            int size_fourier = nb_fields * n_imag * (nlatsLegReduced_ - nlat0_[jm]);
            if (size_fourier > 0) {
                scl_fourier_size = std::max<size_t>(scl_fourier_size, size_fourier);
                scalar_sym_size  = std::max<size_t>(scalar_sym_size, n_imag * nb_fields * size_sym);
                scalar_asym_size = std::max<size_t>(scalar_asym_size, n_imag * nb_fields * size_asym);
            }
        }
    }
    scl_fourier_size = pluto::align_up(scl_fourier_size, alignment/sizeof(double));
    scalar_sym_size  = pluto::align_up(scalar_sym_size,  alignment/sizeof(double));
    scalar_asym_size = pluto::align_up(scalar_asym_size, alignment/sizeof(double));
    scalar_sym_size  = std::max(scalar_sym_size,scalar_asym_size);
    scalar_asym_size = scalar_sym_size;
    double* scalar_sym_buffer;
    double* scalar_asym_buffer;
    double* scl_fourier_sym_buffer;
    double* scl_fourier_asym_buffer;
    {
        ATLAS_TRACE("allocations");
        auto num_threads = static_cast<size_t>(atlas_omp_get_max_threads());
        alloc_aligned(scalar_sym_buffer,       num_threads * scalar_sym_size, "scalar_sym_buffer");
        alloc_aligned(scalar_asym_buffer,      num_threads * scalar_asym_size, "scalar_asym_buffer");
        alloc_aligned(scl_fourier_sym_buffer,  num_threads * scl_fourier_size, "scl_fourier_sym_buffer");
        alloc_aligned(scl_fourier_asym_buffer, num_threads * scl_fourier_size, "scl_fourier_asym_buffer");
    }

    {
        ATLAS_TRACE("Inverse Legendre Transform (GEMM)");
        atlas_omp_parallel_for (int jm = 0; jm <= truncation_; jm++) {
            const int n_imag = (jm ? 2 : 1);
            int size_fourier = nb_fields * n_imag * (nlatsLegReduced_ - nlat0_[jm]);
            if (size_fourier == 0 || jm >= truncation) {
                for (int jfld = 0; jfld < nb_fields; jfld++) {
                    for (int jlat = 0; jlat < nlats; jlat++) {
                        for (int imag = 0; imag < n_imag; imag++) {
                            fourier_out(jfld, jlat, jm, imag) = 0.;
                        }
                    }
                }
                continue;
            }
            else {
                size_t tid{static_cast<size_t>(atlas_omp_get_thread_num())};
                double* scalar_sym       = assume_aligned<alignment>(scalar_sym_buffer       + tid * scalar_sym_size);
                double* scalar_asym      = assume_aligned<alignment>(scalar_asym_buffer      + tid * scalar_asym_size);
                double* scl_fourier_sym  = assume_aligned<alignment>(scl_fourier_sym_buffer  + tid * scl_fourier_size);
                double* scl_fourier_asym = assume_aligned<alignment>(scl_fourier_asym_buffer + tid * scl_fourier_size);
                size_t size_sym  = num_n(truncation_ + 1, jm, true);
                size_t size_asym = num_n(truncation_ + 1, jm, false);

                {
                    ATLAS_TRACE( "Legendre split sym/asym" );
                    idx_t is = 0, ia = 0;
                    const idx_t ioff = (2 * truncation + 3 - jm) * jm / 2 * nb_fields * 2;
                    constexpr auto is_even = [](int n) {
                        // taking bitwise and of n with 1
                        return ((n & 1) == 0);
                    };

                    // the choice between the following two code lines determines whether
                    // total wavenumbers are summed in an ascending or descending order.
                    // The trans library in IFS uses descending order because it should
                    // be more accurate (higher wavenumbers have smaller contributions).
                    // This also needs to be changed when splitting the spectral data in
                    // compute_legendre_polynomials!
                    //for (int jn = truncation+1; jn <= truncation_ + 1; jn++) {
                    for (int jn = truncation_ + 1; jn > truncation; jn--) {
                        if (is_even(jn - jm)) {
                            for (int imag = 0; imag < n_imag; imag++) {
                                for (int jfld = 0; jfld < nb_fields; jfld++) {
                                    scalar_sym[is++] = 0.;
                                }
                            }
                        }
                        else {
                            for (int imag = 0; imag < n_imag; imag++) {
                                for (int jfld = 0; jfld < nb_fields; jfld++) {
                                    scalar_asym[ia++] = 0.;
                                }
                            }
                        }
                    }

                    // the choice between the following two code lines determines whether
                    // total wavenumbers are summed in an ascending or descending order.
                    // The trans library in IFS uses descending order because it should
                    // be more accurate (higher wavenumbers have smaller contributions).
                    // This also needs to be changed when splitting the spectral data in
                    // compute_legendre_polynomials!
                    // for ( int jn = jm; jn <= truncation; jn++ ) {
                    for (int jn = truncation; jn >= jm; jn--) {
                        if (is_even(jn - jm)) {
                            for (int imag = 0; imag < n_imag; imag++) {
                                for (int jfld = 0; jfld < nb_fields; jfld++) {
                                    idx_t idx = jfld + nb_fields * (imag + 2 * (jn - jm));
                                    scalar_sym[is++] = scalar_spectra[idx + ioff];
                                }
                            }
                        }
                        else {
                            for (int imag = 0; imag < n_imag; imag++) {
                                for (int jfld = 0; jfld < nb_fields; jfld++) {
                                    idx_t idx = jfld + nb_fields * (imag + 2 * (jn - jm));
                                    scalar_asym[ia++] = scalar_spectra[idx + ioff];
                                }
                            }
                        }
                    }
                    ATLAS_ASSERT(size_t(ia) == n_imag * nb_fields * size_asym &&
                                 size_t(is) == n_imag * nb_fields * size_sym);
                }
                if (nlatsLegReduced_ - nlat0_[jm] > 0) {
                    ATLAS_TRACE("matrix_multiply (" + std::string(linalg_backend) + ")");
                    {
                        linalg::Matrix A(scalar_sym, nb_fields * n_imag, size_sym);
                        linalg::Matrix B(legendre_sym_ + legendre_sym_begin_[jm] + nlat0_[jm] * size_sym, size_sym,
                                            nlatsLegReduced_ - nlat0_[jm]);
                        linalg::Matrix C(scl_fourier_sym, nb_fields * n_imag, nlatsLegReduced_ - nlat0_[jm]);
                        linalg::matrix_multiply(A, B, C, linalg_backend);
                        /*Log::info() << "sym: ";
                        for ( int j = 0; j < size_sym * ( nlatsLegReduced_ - nlat0_[jm] ); j++ ) {
                            Log::info() << legendre_sym_[j + legendre_sym_begin_[jm] + nlat0_[jm] * size_sym] << " ";
                        }
                        Log::info() << std::endl;*/
                    }
                    if (size_asym > 0) {
                        linalg::Matrix A(scalar_asym, nb_fields * n_imag, size_asym);
                        linalg::Matrix B(legendre_asym_ + legendre_asym_begin_[jm] + nlat0_[jm] * size_asym, size_asym,
                                            nlatsLegReduced_ - nlat0_[jm]);
                        linalg::Matrix C(scl_fourier_asym, nb_fields * n_imag, nlatsLegReduced_ - nlat0_[jm]);
                        linalg::matrix_multiply(A, B, C, linalg_backend);
                        /*Log::info() << "asym: ";
                        for ( int j = 0; j < size_asym * ( nlatsLegReduced_ - nlat0_[jm] ); j++ ) {
                            Log::info() << legendre_asym_[j + legendre_asym_begin_[jm] + nlat0_[jm] * size_asym] << " ";
                        }
                        Log::info() << std::endl;*/
                    }
                    else {
                        for(int i=0; i<nb_fields * n_imag * (nlatsLegReduced_ - nlat0_[jm]); ++i) {
                            scl_fourier_asym[i] = 0.;
                        }
                    }
                }
                {
                    auto posFourier = [&](int jfld, int imag, int jlat, int jm, int nlatsH) {
                        return jfld + nb_fields * (imag + n_imag * (nlatsLegReduced_ - nlat0_[jm] - nlatsH + jlat));
                    };

                    if (nlatsNH_ == nlatsSH_) { // faster
                        int nlatsH = nlatsNH_;

                        ATLAS_TRACE( "Merge symmetric hemispheres" );
                        if (jm == 0) {
                            constexpr int REAL = 0;
                            for (int jlat = 0; jlat < nlatsH; jlat++) {
                                int jlatN = jlat;
                                int jlatS = nlats - jlat - 1;
                                if (nlatsLegReduced_ - nlat0_[jm] - nlatsH + jlat >= 0) {
                                    int idx_jlat = posFourier(0, 0, jlat, jm, nlatsH);
                                    for (int jfld = 0; jfld < nb_fields; jfld++) {
                                        int idx = idx_jlat + jfld;
                                        auto sym  = scl_fourier_sym[idx];
                                        auto asym = scl_fourier_asym[idx];
                                        fourier_out(jfld, jlatN, jm, REAL) = sym + asym;
                                        fourier_out(jfld, jlatS, jm, REAL) = sym - asym;
                                    }
                                }
                                else {
                                    for (int jfld = 0; jfld < nb_fields; jfld++) {
                                        fourier_out(jfld, jlatN, jm, REAL) = 0.;
                                    }
                                }
                            }
                        }
                        else {
                            for (int jlat = 0; jlat < nlatsH; jlat++) {
                                int jlatN = jlat;
                                int jlatS = nlats - jlat - 1;
                                if (nlatsLegReduced_ - nlat0_[jm] - nlatsH + jlat >= 0) {
                                    int idx_jlat = posFourier(0, 0, jlat, jm, nlatsH);
                                    for (int imag = 0; imag < 2; imag++) {
                                        int idx_imag = idx_jlat + imag * nb_fields;
                                        for (int jfld = 0; jfld < nb_fields; jfld++) {
                                            int idx = idx_imag + jfld;
                                            auto sym  = scl_fourier_sym[idx];
                                            auto asym = scl_fourier_asym[idx];
                                            fourier_out(jfld, jlatN, jm, imag) = sym + asym;
                                            fourier_out(jfld, jlatS, jm, imag) = sym - asym;
                                        }
                                    }
                                }
                                else {
                                    for (int jfld = 0; jfld < nb_fields; jfld++) {
                                        for (int imag = 0; imag < 2; imag++) {
                                            fourier_out(jfld, jlatN, jm, imag) = 0.;
                                            fourier_out(jfld, jlatS, jm, imag) = 0.;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else {
                        ATLAS_TRACE( "Merge asymmetric hemispheres" );
                        // northern hemisphere:
                        for (int jlat = 0; jlat < nlatsNH_; jlat++) {
                            if (nlatsLegReduced_ - nlat0_[jm] - nlatsNH_ + jlat >= 0) {
                                int idx_jlat = posFourier(0, 0, jlat, jm, nlatsNH_);
                                for (int imag = 0; imag < n_imag; imag++) {
                                    int idx_imag = idx_jlat + imag * nb_fields;
                                    for (int jfld = 0; jfld < nb_fields; jfld++) {
                                        int idx = idx_imag + jfld;
                                        fourier_out(jfld, jlat, jm, imag) =
                                            scl_fourier_sym[idx] + scl_fourier_asym[idx];
                                    }
                                }
                            }
                            else {
                                for (int jfld = 0; jfld < nb_fields; jfld++) {
                                    for (int imag = 0; imag < n_imag; imag++) {
                                        fourier_out(jfld, jlat, jm, imag) = 0.;
                                    }
                                }
                            }
                            /*for ( int imag = 0; imag < n_imag; imag++ ) {
                            for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                                if ( scl_fourier[posMethod( jfld, imag, jlat, jm, nb_fields, nlats )] > 0. ) {
                                    Log::info() << "jm=" << jm << " jlat=" << jlat << " nlatsLeg_=" << nlatsLeg_
                                                << " nlat0=" << nlat0_[jm] << " nlatsNH=" << nlatsNH_ << std::endl;
                                }
                            }*/
                        }
                        // southern hemisphere:
                        for (int jlat = 0; jlat < nlatsSH_; jlat++) {
                            int jslat = nlats - jlat - 1;
                            if (nlatsLegReduced_ - nlat0_[jm] - nlatsSH_ + jlat >= 0) {
                                int idx_jlat = posFourier(0, 0, jlat, jm, nlatsSH_);
                                for (int imag = 0; imag < n_imag; imag++) {
                                    int idx_imag = idx_jlat + imag * nb_fields;
                                    for (int jfld = 0; jfld < nb_fields; jfld++) {
                                        int idx = idx_imag + jfld;
                                        fourier_out(jfld, jslat, jm, imag) =
                                            scl_fourier_sym[idx] - scl_fourier_asym[idx];
                                    }
                                }
                            }
                            else {
                                for (int jfld = 0; jfld < nb_fields; jfld++) {
                                    for (int imag = 0; imag < n_imag; imag++) {
                                        fourier_out(jfld, jslat, jm, imag) = 0.;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    {
        ATLAS_TRACE("deallocations");
        auto num_threads = static_cast<size_t>(atlas_omp_get_max_threads());
        free_aligned(scalar_sym_buffer,       num_threads * scalar_sym_size,  "scalar_sym_buffer");
        free_aligned(scalar_asym_buffer,      num_threads * scalar_asym_size, "scalar_asym_buffer");
        free_aligned(scl_fourier_sym_buffer,  num_threads * scl_fourier_size, "scl_fourier_sym_buffer");
        free_aligned(scl_fourier_asym_buffer, num_threads * scl_fourier_size, "scl_fourier_asym_buffer");
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_fourier_regular(const int nlats, const int nlons, const int nb_fields, double scl_fourier[],
                                          double gp_fields[], const eckit::Configuration&) const {
    // Fourier transformation:
    if (useFFT_) {
        int num_complex = (nlonsMaxGlobal_ / 2) + 1;
        {
            ATLAS_TRACE("Inverse Fourier Transform ("+fft_->type()+", RegularGrid)");
            for (int jfld = 0; jfld < nb_fields; jfld++) {
                atlas_omp_parallel_for (int jlat = 0; jlat < nlats; jlat++) {
                    constexpr int REAL = 0;
                    constexpr int IMAG = 1;
                    auto fft_data_in = fft_data_->one_in(jlat);
                    int m_trunc = std::min(num_complex,truncation_+1);
                    if (m_trunc > 0) {
                        int jm = 0;
                        fft_data_in[jm] = std::complex<double>{
                            scl_fourier[posMethod(jfld, REAL, jlat, jm, nb_fields, nlats)],
                            0.
                        };
                    }
                    for (int jm= 1; jm < m_trunc; jm++) {
                        fft_data_in[jm] = std::complex<double>{
                            scl_fourier[posMethod(jfld, REAL, jlat, jm, nb_fields, nlats)],
                            scl_fourier[posMethod(jfld, IMAG, jlat, jm, nb_fields, nlats)]
                        };
                    }
                    for(int jm = m_trunc; jm < num_complex; ++jm) {
                        fft_data_in[jm] = std::complex<double>{0., 0.};
                    }
                }

                fft_->inverse_c2r_many( fft_data_->many_in(), fft_data_->many_out() );

                atlas_omp_parallel_for (int jlat = 0; jlat < nlats; jlat++) {
                    auto fft_data_out = fft_data_->one_out(jlat);
                    for (int jlon = 0; jlon < nlons; jlon++) {
                        int j = jlon + jlonMin_[0];
                        if (j >= nlonsMaxGlobal_) {
                            j -= nlonsMaxGlobal_;
                        }
                        gp_fields[jlon + nlons * (jlat + nlats * jfld)] = fft_data_out[j];
                    }
                }
            }
        }
    }
    else {
        linalg::dense::Backend linalg_backend{linalg_backend_};
        // dgemm-method 1
        {
            ATLAS_TRACE("Inverse Fourier Transform (NoFFT,matrix_multiply=" + detect_linalg_backend(linalg_backend_) +
                        ")");
            linalg::Matrix A(fourier_, nlons, (truncation_ + 1) * 2);
            linalg::Matrix B(scl_fourier, (truncation_ + 1) * 2, nb_fields * nlats);
            linalg::Matrix C(gp_fields, nlons, nb_fields * nlats);

            // TODO: OpenMP,  Either use OpenMP backend or we manually multithread?
            linalg::matrix_multiply(A, B, C, linalg_backend);
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_fourier_reduced(const int nlats, const StructuredGrid& g, const int nb_fields,
                                          double scl_fourier[], double gp_fields[], const eckit::Configuration&) const {
    // Fourier transformation:
    if (useFFT_) {
        ATLAS_TRACE("Inverse Fourier Transform ("+fft_->type()+", ReducedGrid)");
        std::vector<size_t, pluto::allocator<size_t>> gp_offsets(nlats+1, memory_resource());
        for (int jlat = 0; jlat < nlats; jlat++) {
            gp_offsets[jlat + 1] = gp_offsets[jlat] + g.nx(jlat);
        }
        constexpr int scl_fourier_stride_jm = 2;
        const int scl_fourier_stride_jlat = (truncation_ + 1) * scl_fourier_stride_jm;
        const int scl_fourier_stride_jfld = nlats * scl_fourier_stride_jlat;
        for (int jfld = 0; jfld < nb_fields; jfld++) {
            const int scl_fourier_offset_jfld = scl_fourier_stride_jfld * jfld;
            auto scl_fourier_index = [scl_fourier_offset_jfld, scl_fourier_stride_jlat](int jlat, int jm) {
                return scl_fourier_offset_jfld
                     + scl_fourier_stride_jlat * jlat 
                     + scl_fourier_stride_jm * jm;
            };
            atlas_omp_parallel_for (int jlat = 0; jlat < nlats; jlat++) {
                int tid = atlas_omp_get_thread_num();
                auto fft_data_in  = assume_aligned<alignment>(fft_data_->one_in(tid));
                auto fft_data_out = assume_aligned<alignment>(fft_data_->one_out(tid));

                // Load fft_data_in
                {
                    constexpr int REAL = 0;
                    constexpr int IMAG = 1;
                    int num_complex     = (nlonsGlobal_[jlat] / 2) + 1;
                    int m_trunc = std::min(num_complex,truncation_+1);
                    if (m_trunc > 0) {
                        int jm = 0;
                        fft_data_in[jm] = std::complex<double>{
                            scl_fourier[scl_fourier_index(jlat,jm)],
                            0.
                        };
                    }
                    for (int jm = 1; jm < m_trunc; jm++) {
                        int idx = scl_fourier_index(jlat,jm);
                        fft_data_in[jm] = std::complex<double>{
                            scl_fourier[idx+REAL],
                            scl_fourier[idx+IMAG]
                        };
                    }
                    for (int jm = m_trunc; jm < num_complex; jm++) {
                        fft_data_in[jm] = std::complex<double>{0., 0.};
                    }
                }

                fft_->inverse_c2r(nlonsGlobal_[jlat], fft_data_in, fft_data_out);

                // Unload fft_data_out
                {
                    int jgp = gp_offsets[jlat] + jfld * gp_offsets.back();

                    if (jlonMin_[jlat] == 0) {
                        for (int jlon = 0; jlon < g.nx(jlat); jlon++) {
                            gp_fields[jgp++] = fft_data_out[jlon];
                        }
                    }
                    else {
                        for (int jlon = 0; jlon < g.nx(jlat); jlon++) {
                            int j = jlon + jlonMin_[jlat];
                            if (j >= nlonsGlobal_[jlat]) {
                                j -= nlonsGlobal_[jlat];
                            }
                            gp_fields[jgp++] = fft_data_out[j];
                        }
                    }
                }
            }
        }
    }
    else {
        throw_NotImplemented(
            "Using dgemm in Fourier transform for reduced grids is extremely slow.",
            Here());
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_unstructured_precomp(const int truncation, const int nb_fields, const int nb_vordiv_fields,
                                               const double scalar_spectra[], double gp_fields[],
                                               const eckit::Configuration&) const {
    ATLAS_TRACE("invtrans_uv unstructured");

    const int nlats        = grid_.size();
    const int size_fourier = nb_fields * 2;
    double* scl_fourier;
    double* scl_fourier_tp;
    double* fouriertp;
    double* gp_opt;
    size_t scl_fourier_size = size_fourier * (truncation)*nlats;
    size_t scl_fourier_tp_size =  size_fourier * (truncation);
    size_t fouriertp_size = 2 * (truncation);
    size_t gp_opt_size = nb_fields;
    alloc_aligned(scl_fourier, scl_fourier_size, "scl_fourier");
    alloc_aligned(scl_fourier_tp, scl_fourier_tp_size, "scl_fourier_tp");
    alloc_aligned(fouriertp, fouriertp_size, "fouriertp");
    alloc_aligned(gp_opt, gp_opt_size, "gp_opt");

    linalg::dense::Backend linalg_backend{linalg_backend_};

    {
        ATLAS_TRACE("Inverse Legendre Transform (GEMM," + std::string(linalg_backend) + ")");
        atlas_omp_parallel_for (int jm = 0; jm < truncation; jm++) {
            const int noff = (2 * truncation + 3 - jm) * jm / 2, ns = truncation - jm + 1;
            linalg::Matrix A(
                eckit::linalg::Matrix(const_cast<double*>(scalar_spectra) + nb_fields * 2 * noff, nb_fields * 2, ns));
            linalg::Matrix B(legendre_ + noff * nlats, ns, nlats);
            linalg::Matrix C(scl_fourier + jm * size_fourier * nlats, nb_fields * 2, nlats);
            linalg::matrix_multiply(A, B, C, linalg_backend);
        }
    }

    // loop over all points:
    {
        ATLAS_TRACE("Inverse Fourier Transform (NoFFT," + std::string(linalg_backend) + ")");
        int ip = 0;
// TODO: OpenMP
        for (const PointLonLat p : grid_.lonlat()) {
            const double lon = p.lon() * util::Constants::degreesToRadians();
            const double lat = p.lat() * util::Constants::degreesToRadians();
            {
                //ATLAS_TRACE( "opt transposition in Fourier" );
                for (int jm = 0; jm < truncation; jm++) {
                    int idx = nb_fields * 2 * (ip + nlats * jm);
                    for (int imag = 0; imag < 2; imag++) {
                        for (int jfld = 0; jfld < nb_fields; jfld++) {
                            const int pos_tp = imag + 2 * (jm + (truncation) * (jfld));
                            //int pos  = jfld + nb_fields * ( imag + 2 * ( jm ) );
                            scl_fourier_tp[pos_tp] = scl_fourier[idx++];  // = scl_fourier[pos]
                        }
                    }
                }
            }

            // Fourier transformation:
            {
                //ATLAS_TRACE( "opt compute fouriertp" );
                int idx          = 0;
                fouriertp[idx++] = 1.;  // real part
                fouriertp[idx++] = 0.;  // imaginary part
                for (int jm = 1; jm < truncation; jm++) {
                    fouriertp[idx++] = +2. * std::cos(jm * lon);  // real part
                    fouriertp[idx++] = -2. * std::sin(jm * lon);  // imaginary part
                }
            }
            {
                //ATLAS_TRACE( "opt Fourier dgemm" );
                linalg::Matrix A(fouriertp, 1, (truncation)*2);
                linalg::Matrix B(scl_fourier_tp, (truncation)*2, nb_fields);
                linalg::Matrix C(gp_opt, 1, nb_fields);
                linalg::matrix_multiply(A, B, C, linalg_backend);
                for (int j = 0; j < nb_fields; j++) {
                    gp_fields[ip + j * grid_.size()] = gp_opt[j];
                }
            }
            // Computing u,v from U,V:
            {
                if (nb_vordiv_fields > 0) {
                    //ATLAS_TRACE( " u,v from U,V" );
                    double coslat = std::cos(lat);
                    for (int j = 0; j < 2 * nb_vordiv_fields && j < nb_fields; j++) {
                        gp_fields[ip + j * grid_.size()] /= coslat;
                    }
                }
            }
            ++ip;
        }
    }
    free_aligned(scl_fourier, scl_fourier_size, "scl_fourier");
    free_aligned(scl_fourier_tp, scl_fourier_tp_size, "scl_fourier");
    free_aligned(fouriertp, fouriertp_size, "scl_fourier");
    free_aligned(gp_opt, gp_opt_size, "scl_fourier");
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_unstructured(const int truncation, const int nb_fields, const int nb_vordiv_fields,
                                       const double scalar_spectra[], double gp_fields[],
                                       const eckit::Configuration& config) const {
    ATLAS_TRACE("invtrans_unstructured");

    if (warning(config)) {
        Log::warning() << "WARNING: Spectral transforms could take a long time (unstructured grid approach). Results "
                          "may contain aliasing errors."
                       << std::endl;
    }

    linalg::dense::Backend linalg_backend{linalg_backend_};

    double* zfn;
    size_t zfn_size = (truncation + 1) * (truncation + 1);
    alloc_aligned(zfn, zfn_size, "zfn");
    compute_zfn(truncation, zfn);
    int size_fourier = nb_fields * 2;
    double* legendre;
    double* scl_fourier;
    double* scl_fourier_tp;
    double* fouriertp;
    double* gp_opt;

    size_t _legendre_size = legendre_size(truncation + 1);
    size_t scl_fourier_size = size_fourier * (truncation + 1);
    size_t scl_fourier_tp_size =  size_fourier * (truncation + 1);
    size_t fouriertp_size = 2 * (truncation + 1);
    size_t gp_opt_size = nb_fields;

    alloc_aligned(legendre, legendre_size(truncation + 1), "legendre");
    alloc_aligned(scl_fourier, scl_fourier_size, "scl_fourier");
    alloc_aligned(scl_fourier_tp, scl_fourier_tp_size, "scl_fourier_tp");
    alloc_aligned(fouriertp, fouriertp_size, "fouriertp");
    alloc_aligned(gp_opt, gp_opt_size, "gp_opt");


    // loop over all points:
    int ip = 0;
    LegendrePolynomialsWorkspace w{truncation};
    // TODO OpenMP
    for (const PointLonLat p : grid_.lonlat()) {
        const double lon = p.lon() * util::Constants::degreesToRadians();
        const double lat = p.lat() * util::Constants::degreesToRadians();
        compute_legendre_polynomials_lat(truncation, lat, legendre, zfn, w);
        // Legendre transform:
        {
            //ATLAS_TRACE( "opt Legendre dgemm" );
            ATLAS_TRACE("Legendre matrix_multiply (" + std::string(linalg_backend) + ")");
            for (int jm = 0; jm <= truncation; jm++) {
                const int noff = (2 * truncation + 3 - jm) * jm / 2, ns = truncation - jm + 1;
                linalg::Matrix A(eckit::linalg::Matrix(const_cast<double*>(scalar_spectra) + nb_fields * 2 * noff,
                                                       nb_fields * 2, ns));
                linalg::Matrix B(legendre + noff, ns, 1);
                linalg::Matrix C(scl_fourier + jm * size_fourier, nb_fields * 2, 1);
                linalg::matrix_multiply(A, B, C, linalg_backend);
            }
        }
        {
            //ATLAS_TRACE( "opt transposition in Fourier" );
            int idx = 0;
            for (int jm = 0; jm < truncation + 1; jm++) {
                for (int imag = 0; imag < 2; imag++) {
                    for (int jfld = 0; jfld < nb_fields; jfld++) {
                        const int pos_tp = imag + 2 * (jm + (truncation + 1) * (jfld));
                        //int pos  = jfld + nb_fields * ( imag + 2 * ( jm ) );
                        scl_fourier_tp[pos_tp] = scl_fourier[idx++];  // = scl_fourier[pos]
                    }
                }
            }
        }

        // Fourier transformation:
        int idx          = 0;
        fouriertp[idx++] = 1.;  // real part
        fouriertp[idx++] = 0.;  // imaginary part
        for (int jm = 1; jm < truncation + 1; jm++) {
            fouriertp[idx++] = +2. * std::cos(jm * lon);  // real part
            fouriertp[idx++] = -2. * std::sin(jm * lon);  // imaginary part
        }
        {
            ATLAS_TRACE("Fourier matrix_multiply (" + std::string(linalg_backend) + ")");
            linalg::Matrix A(fouriertp, 1, (truncation + 1) * 2);
            linalg::Matrix B(scl_fourier_tp, (truncation + 1) * 2, nb_fields);
            linalg::Matrix C(gp_opt, 1, nb_fields);
            linalg::matrix_multiply(A, B, C, linalg_backend);
            for (int j = 0; j < nb_fields; j++) {
                gp_fields[ip + j * grid_.size()] = gp_opt[j];
            }
        }
        // Computing u,v from U,V:
        {
            if (nb_vordiv_fields > 0) {
                //ATLAS_TRACE( "u,v from U,V" );
                const double coslat = std::cos(lat);
                for (int j = 0; j < 2 * nb_vordiv_fields && j < nb_fields; j++) {
                    gp_fields[ip + j * grid_.size()] /= coslat;
                }
            }
        }
        ++ip;
    }
    free_aligned(legendre, _legendre_size, "legendre");
    free_aligned(scl_fourier, scl_fourier_size, "scl_fourier");
    free_aligned(scl_fourier_tp, scl_fourier_tp_size, "scl_fourier_tp");
    free_aligned(fouriertp, fouriertp_size, "fouriertp");
    free_aligned(gp_opt, gp_opt_size, "gp_opt");
    free_aligned(zfn, zfn_size, "zfn");
}

//-----------------------------------------------------------------------------
// Routine to compute the spectral transform by using a Local Fourier transformation
// for a grid (same latitude for all longitudes, allows to compute Legendre functions
// once for all longitudes). U and v components are divided by cos(latitude) for
// nb_vordiv_fields > 0.
//
// Legendre polynomials are computed up to truncation_+1 to be accurate for vorticity and
// divergence computation. The parameter truncation is the truncation used in storing the
// spectral data scalar_spectra and can be different from truncation_. If truncation is
// larger than truncation_+1 the transform will behave as if the spectral data was truncated
// to truncation_+1.
//
// Author:
// Andreas Mueller *ECMWF*
//
void TransLocal::invtrans_uv(const int truncation, const int nb_scalar_fields, const int nb_vordiv_fields,
                             const double scalar_spectra[], double gp_fields[],
                             const eckit::Configuration& config) const {
    if (nb_scalar_fields > 0) {
        int nb_fields = nb_scalar_fields;

        // Transform
        if (StructuredGrid(grid_) && not grid_.projection()) {
            auto g = StructuredGrid(grid_);
            ATLAS_TRACE("invtrans_uv structured");
            int nlats            = g.ny();
            int nlons            = g.nxmax();
            double* scl_fourier;
            size_t scl_fourier_size = size_t(nb_fields) * size_t(nlats) * size_t(truncation_ + 1) * 2;
            alloc_aligned(scl_fourier, scl_fourier_size, "scl_fourier");

            // ATLAS-159 workaround begin
            atlas_omp_parallel_for (int i = 0; i < scl_fourier_size; ++i) {
                scl_fourier[i] = 0.;
            }
            // ATLAS-159 workaround end

            // Legendre transformation:
            invtrans_legendre(truncation, nlats, nb_scalar_fields, nb_vordiv_fields, scalar_spectra, scl_fourier,
                              config);

            // Fourier transformation:
            if (RegularGrid(gridGlobal_)) {
                invtrans_fourier_regular(nlats, nlons, nb_fields, scl_fourier, gp_fields, config);
            }
            else {
                invtrans_fourier_reduced(nlats, g, nb_fields, scl_fourier, gp_fields, config);
            }

            // Computing u,v from U,V:
            {
                if (nb_vordiv_fields > 0) {
                    ATLAS_TRACE("compute u,v from U,V");
                    std::vector<double, pluto::allocator<double>> coslatinvs(nlats, memory_resource());
                    for (idx_t j = 0; j < nlats; ++j) {
                        double lat = g.y(j);
                        if (lat > latPole) {
                            lat = latPole;
                        }
                        if (lat < -latPole) {
                            lat = -latPole;
                        }
                        double coslat = std::cos(lat * util::Constants::degreesToRadians());
                        coslatinvs[j] = 1. / coslat;
                        //Log::info() << "lat=" << g.y( j ) << " coslat=" << coslat << std::endl;
                    }
                    int idx = 0;
                    for (idx_t jfld = 0; jfld < 2 * nb_vordiv_fields && jfld < nb_fields; jfld++) {
                        for (idx_t jlat = 0; jlat < g.ny(); jlat++) {
                            for (idx_t jlon = 0; jlon < g.nx(jlat); jlon++) {
                                gp_fields[idx] *= coslatinvs[jlat];
                                idx++;
                            }
                        }
                    }
                }
            }
            free_aligned(scl_fourier, scl_fourier_size, "scl_fourier");
        }
        else {
            if (unstruct_precomp_) {
                invtrans_unstructured_precomp(truncation, nb_scalar_fields, nb_vordiv_fields, scalar_spectra, gp_fields,
                                              config);
            }
            else {
                invtrans_unstructured(truncation, nb_scalar_fields, nb_vordiv_fields, scalar_spectra, gp_fields,
                                      config);
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans(const int nb_vordiv_fields, const double vorticity_spectra[],
                          const double divergence_spectra[], double gp_fields[],
                          const eckit::Configuration& config) const {
    invtrans(0, nullptr, nb_vordiv_fields, vorticity_spectra, divergence_spectra, gp_fields, config);
}

// --------------------------------------------------------------------------------------------------------------------

void extend_truncation(const int old_truncation, const int nb_fields, const double old_spectra[], double new_spectra[]) {

    int new_truncation = old_truncation + 1;
    int old_size = 2 * legendre_size(old_truncation) * nb_fields;
    int new_size = 2 * legendre_size(new_truncation) * nb_fields;

    int k = 0, k_old = 0;
    for (int m = 0; m <= new_truncation; m++) {             // zonal wavenumber
        for (int n = m; n <= new_truncation; n++) {         // total wavenumber
            for (int imag = 0; imag < 2; imag++) {              // imaginary/real part
                for (int jfld = 0; jfld < nb_fields; jfld++) {  // field
                    if (m == new_truncation || n == new_truncation) {
                        new_spectra[k++] = 0.;
                    }
                    else {
                        new_spectra[k++] = old_spectra[k_old++];
                    }
                }
            }
        }
    }
    ATLAS_ASSERT(k==new_size);
    ATLAS_ASSERT(k_old==old_size);
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans(const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                          const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                          const eckit::Configuration& config) const {
    int nb_gp = grid_.size();
    if (nb_vordiv_fields > 0) {
        // collect all spectral data into one array "all_spectra":
        ATLAS_TRACE("TransLocal::invtrans");
        int nb_vordiv_spec_ext = 2 * legendre_size(truncation_ + 1) * nb_vordiv_fields;
        std::vector<double,pluto::allocator<double>> U_ext(memory_resource());
        std::vector<double,pluto::allocator<double>> V_ext(memory_resource());
        std::vector<double,pluto::allocator<double>> scalar_ext(memory_resource());
        if (nb_vordiv_fields > 0) {
            std::vector<double,pluto::allocator<double>> vorticity_spectra_extended(nb_vordiv_spec_ext, memory_resource());
            std::vector<double,pluto::allocator<double>> divergence_spectra_extended(nb_vordiv_spec_ext, memory_resource());
            U_ext.resize(nb_vordiv_spec_ext);
            V_ext.resize(nb_vordiv_spec_ext);

            {
                ATLAS_TRACE("extend vordiv");
                // increase truncation in vorticity_spectra and divergence_spectra:
                extend_truncation(truncation_, nb_vordiv_fields, vorticity_spectra, vorticity_spectra_extended.data());
                extend_truncation(truncation_, nb_vordiv_fields, divergence_spectra,
                                  divergence_spectra_extended.data());
            }

            {
                ATLAS_TRACE("vordiv to UV");
                // call vd2uv to compute u and v in spectral space
                trans::VorDivToUV vordiv_to_UV_ext(truncation_ + 1, option::type("local"));
                vordiv_to_UV_ext.execute(nb_vordiv_spec_ext, nb_vordiv_fields, vorticity_spectra_extended.data(),
                                         divergence_spectra_extended.data(), U_ext.data(), V_ext.data());
            }
        }
        if (nb_scalar_fields > 0) {
            int nb_scalar_ext = 2 * legendre_size(truncation_ + 1) * nb_scalar_fields;
            scalar_ext.resize(nb_scalar_ext);
            extend_truncation(truncation_, nb_scalar_fields, scalar_spectra, scalar_ext.data());
        }
        int nb_all_fields = 2 * nb_vordiv_fields + nb_scalar_fields;
        int nb_all_size   = 2 * legendre_size(truncation_ + 1) * nb_all_fields;
        std::vector<double, pluto::allocator<double>> all_spectra(nb_all_size, memory_resource());
        int k = 0, i = 0, j = 0, l = 0;
        {
            ATLAS_TRACE("merge all spectra");
            for (int m = 0; m <= truncation_ + 1; m++) {                       // zonal wavenumber
                for (int n = m; n <= truncation_ + 1; n++) {                   // total wavenumber
                    for (int imag = 0; imag < 2; imag++) {                     // imaginary/real part
                        for (int jfld = 0; jfld < nb_vordiv_fields; jfld++) {  // vorticity fields
                            all_spectra[k++] = U_ext[i++];
                        }
                        for (int jfld = 0; jfld < nb_vordiv_fields; jfld++) {  // divergence fields
                            all_spectra[k++] = V_ext[j++];
                        }
                        for (int jfld = 0; jfld < nb_scalar_fields; jfld++) {  // scalar fields
                            all_spectra[k++] = scalar_ext[l++];
                        }
                    }
                }
            }
        }
        int nb_vordiv_size = 2 * legendre_size(truncation_ + 1) * nb_vordiv_fields;
        int nb_scalar_size = 2 * legendre_size(truncation_ + 1) * nb_scalar_fields;
        ATLAS_ASSERT(k == nb_all_size);
        ATLAS_ASSERT(i == nb_vordiv_size);
        ATLAS_ASSERT(j == nb_vordiv_size);
        ATLAS_ASSERT(l == nb_scalar_size);
        invtrans_uv(truncation_ + 1, nb_all_fields, nb_vordiv_fields, all_spectra.data(), gp_fields, config);
    }
    else {
        if (nb_scalar_fields > 0) {
            invtrans_uv(truncation_, nb_scalar_fields, 0, scalar_spectra, gp_fields + 2 * nb_gp * nb_vordiv_fields,
                        config);
        }
    }
}

void TransLocal::invtrans_adj(const int nb_scalar_fields, const double gp_fields[], double scalar_spectra[],
                              const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

void TransLocal::invtrans_adj(const int nb_vordiv_fields, const double gp_fields[], double vorticity_spectra[],
                              double divergence_spectra[], const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

void TransLocal::invtrans_adj(const int nb_scalar_fields, const double gp_fields[], const int nb_vordiv_fields,
                              double vorticity_spectra[], double divergence_spectra[], double scalar_spectra[],
                              const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans(const Field& gpfield, Field& spfield, const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans(const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans_wind2vordiv(const Field& gpwind, Field& spvor, Field& spdiv,
                                      const eckit::Configuration& config) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}


void TransLocal::dirtrans_adj(const Field& spfield, Field& gpfield,
                              const eckit::Configuration& config) const {
  ATLAS_NOTIMPLEMENTED;
  // Not implemented and not planned.
  // Use the TransIFS implementation instead.
}

void TransLocal::dirtrans_adj(const FieldSet& spfields, FieldSet& gpfields,
                              const eckit::Configuration& config) const {
  ATLAS_NOTIMPLEMENTED;
  // Not implemented and not planned.
  // Use the TransIFS implementation instead.
}

void TransLocal::dirtrans_wind2vordiv_adj(const Field& spvor, const Field& spdiv, Field& gpwind,
                                          const eckit::Configuration& config) const {
  ATLAS_NOTIMPLEMENTED;
  // Not implemented and not planned.
  // Use the TransIFS implementation instead.
}


// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans(const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                          const eckit::Configuration&) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans(const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                          double divergence_spectra[], const eckit::Configuration&) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
