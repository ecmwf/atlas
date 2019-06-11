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
#include "eckit/linalg/Matrix.h"
#include "eckit/log/Bytes.h"
#include "eckit/parser/JSON.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/VorDivToUV.h"
#include "atlas/trans/detail/TransFactory.h"
#include "atlas/trans/local/LegendrePolynomials.h"
#include "atlas/util/Constants.h"

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
static TransBuilderGrid<TransLocal> builder( "local", "local" );
}  // namespace

namespace {
class TransParameters {
public:
    TransParameters( const eckit::Configuration& config ) : config_( config ) {}
    ~TransParameters() {}

    /*
     * For the future
     */
    //    bool scalar_derivatives() const { return config_.getBool( "scalar_derivatives", false ); }

    //    bool wind_EW_derivatives() const { return config_.getBool( "wind_EW_derivatives", false ); }

    //    bool vorticity_divergence_fields() const { return config_.getBool( "vorticity_divergence_fields", false ); }

    //    std::string read_legendre() const { return config_.getString( "read_legendre", "" ); }
    //    bool global() const { return config_.getBool( "global", false ); }

    //    std::string read_fft() const { return config_.getString( "read_fft", "" ); }


    std::string write_legendre() const { return config_.getString( "write_legendre", "" ); }

    std::string write_fft() const { return config_.getString( "write_fft", "" ); }


    bool export_legendre() const { return config_.getBool( "export_legendre", false ); }

    int warning() const { return config_.getInt( "warning", 1 ); }

    int fft() const {
        static const std::map<std::string, int> string_to_FFT = {{"OFF", static_cast<int>( option::FFT::OFF )},
                                                                 {"FFTW", static_cast<int>( option::FFT::FFTW )}};
#ifdef ATLAS_HAVE_FFTW
        std::string fft_default = "FFTW";
#else
        std::string fft_default = "OFF";
#endif
        return string_to_FFT.at( config_.getString( "fft", fft_default ) );
    }

private:
    const eckit::Configuration& config_;
};

struct ReadCache {
    ReadCache( const void* cache ) {
        begin = reinterpret_cast<const char*>( cache );
        pos   = 0;
    }
    template <typename T>
    T* read( size_t size ) {
        const T* v = reinterpret_cast<const T*>( begin + pos );
        pos += size * sizeof( T );
        return const_cast<T*>( v );
    }

    const char* begin;
    size_t pos;
};

struct WriteCache {
    WriteCache( const eckit::PathName& file_path ) : dh_( file_path.fileHandle( /*overwrite = */ true ) ) {
        if ( file_path.exists() ) {
            std::stringstream err;
            err << "Cannot open cache file " << file_path << " for writing as it already exists. Remove first.";
            throw_Exception( err.str(), Here() );
        }
        dh_->openForWrite( 0 );
        pos = 0;
    }
    ~WriteCache() { dh_->close(); }
    template <typename T>
    void write( const T* v, long size ) {
        dh_->write( v, size * sizeof( T ) );
        pos += size * sizeof( T );
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

size_t legendre_size( const size_t truncation ) {
    return ( truncation + 2 ) * ( truncation + 1 ) / 2;
}

//int nlats_northernHemisphere( const int nlats ) {
//    return ceil( nlats / 2. );
//    // using ceil here should make it possible to have odd number of latitudes (with the centre latitude being the equator)
//}

size_t num_n( const int truncation, const int m, const bool symmetric ) {
    int len = ( truncation - m + ( symmetric ? 2 : 1 ) ) / 2;
    ATLAS_ASSERT( len >= 0 );
    return size_t( len );
}


[[noreturn]] void throw_AllocationFailed( size_t bytes, const eckit::CodeLocation& loc ) {
    std::stringstream ss;
    ss << "AllocationFailed: Could not allocate " << eckit::Bytes( bytes );
    throw_Exception( ss.str(), loc );
}


void alloc_aligned( double*& ptr, size_t n ) {
    const size_t alignment = 64 * sizeof( double );
    size_t bytes           = sizeof( double ) * n;
    int err                = posix_memalign( (void**)&ptr, alignment, bytes );
    if ( err ) { throw_AllocationFailed( bytes, Here() ); }
}

void free_aligned( double*& ptr ) {
    free( ptr );
    ptr = nullptr;
}

void alloc_aligned( double*& ptr, size_t n, const char* msg ) {
    ATLAS_ASSERT( msg );
    Log::debug() << "TransLocal: allocating '" << msg << "': " << eckit::Bytes( sizeof( double ) * n ) << std::endl;
    alloc_aligned( ptr, n );
}

void free_aligned( double*& ptr, const char* msg ) {
    ATLAS_ASSERT( msg );
    Log::debug() << "TransLocal: dellocating '" << msg << "'" << std::endl;
    free_aligned( ptr );
}

size_t add_padding( size_t n ) {
    return size_t( std::ceil( n / 8. ) ) * 8;
}

}  // namespace

int fourier_truncation( const int truncation,    // truncation
                        const int nx,            // number of longitudes
                        const int /*nxmax*/,     // maximum nx
                        const int ndgl,          // number of latitudes
                        const double lat,        // latitude in radian
                        const bool fullgrid ) {  // regular grid
    int trc     = truncation;
    int trclin  = ndgl - 1;
    int trcquad = ndgl * 2 / 3 - 1;
    if ( truncation >= trclin || fullgrid ) {
        // linear
        trc = ( nx - 1 ) / 2;
    }
    else if ( truncation >= trcquad ) {
        // quadratic
        double weight = 3 * ( trclin - truncation ) / ndgl;
        double sqcos  = std::pow( std::cos( lat ), 2 );

        trc = static_cast<int>( ( nx - 1 ) / ( 2 + weight * sqcos ) );
    }
    else {
        // cubic
        double sqcos = std::pow( std::cos( lat ), 2 );

        trc = static_cast<int>( ( nx - 1 ) / ( 2 + sqcos ) - 1 );
    }
    trc = std::min( truncation, trc );
    return trc;
}

namespace detail {
struct FFTW_Data {
#if ATLAS_HAVE_FFTW
    fftw_complex* in;
    double* out;
    std::vector<fftw_plan> plans;
#endif
};
}  // namespace detail


// --------------------------------------------------------------------------------------------------------------------
// Class TransLocal
// --------------------------------------------------------------------------------------------------------------------

const eckit::linalg::LinearAlgebra& linear_algebra_backend() {
    if ( eckit::linalg::LinearAlgebra::hasBackend( "mkl" ) ) {
        return eckit::linalg::LinearAlgebra::getBackend( "mkl" );
    }
    // Default backend
    return eckit::linalg::LinearAlgebra::backend();
}

bool TransLocal::warning( const eckit::Configuration& config ) const {
    int warning = warning_;
    config.get( "warning", warning );
    return ( warning > 0 && grid_.size() >= warning );
}

TransLocal::TransLocal( const Cache& cache, const Grid& grid, const Domain& domain, const long truncation,
                        const eckit::Configuration& config ) :
    grid_( grid, domain ),
    truncation_( static_cast<int>( truncation ) ),
    precompute_( config.getBool( "precompute", true ) ),
    cache_( cache ),
    legendre_cache_( cache.legendre().data() ),
    legendre_cachesize_( cache.legendre().size() ),
    fft_cache_( cache.fft().data() ),
    fft_cachesize_( cache.fft().size() ),
    fftw_( new detail::FFTW_Data ),
    linalg_( linear_algebra_backend() ),
    warning_( TransParameters( config ).warning() ) {
    ATLAS_TRACE( "TransLocal constructor" );
    double fft_threshold = 0.0;  // fraction of latitudes of the full grid down to which FFT is used.
    // This threshold needs to be adjusted depending on the dgemm and FFT performance of the machine
    // on which this code is running!
    int nlats         = 0;
    int nlonsMax      = 0;
    int neqtr         = 0;
    useFFT_           = TransParameters( config ).fft();
    unstruct_precomp_ = ( config.has( "precompute" ) ? precompute_ : false );
    no_symmetry_      = false;
    nlatsNH_          = 0;
    nlatsSH_          = 0;
    nlatsLeg_         = 0;
    nlatsLegDomain_   = 0;
    nlatsLegReduced_  = 0;
    bool useGlobalLeg = true;
    bool no_nest      = false;

    if ( StructuredGrid( grid_ ) && not grid_.projection() ) {
        StructuredGrid g( grid_ );
        nlats    = g.ny();
        nlonsMax = g.nxmax();

        // check location of domain relative to the equator:
        for ( idx_t j = 0; j < nlats; ++j ) {
            // assumptions: latitudes in g.y(j) are monotone and decreasing
            // no assumption on whether we have 0, 1 or 2 latitudes at the equator
            double lat = g.y( j );
            ( eckit::types::is_approximately_equal( lat, 0. ) ? neqtr : ( lat < 0 ? nlatsSH_ : nlatsNH_ ) )++;
        }
        if ( neqtr > 0 ) {
            nlatsNH_++;
            nlatsSH_++;
        }
        if ( nlatsNH_ >= nlatsSH_ ) { nlatsLegDomain_ = nlatsNH_; }
        else {
            nlatsLegDomain_ = nlatsSH_;
        }

        gridGlobal_ = grid;
        if ( not gridGlobal_.domain().global() ) {
            // The grid is not a nest of a global grid
            if ( RegularGrid( grid_ ) ) {
                no_nest         = true;
                no_symmetry_    = true;
                useFFT_         = false;
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
                throw_NotImplemented( log.str(), Here() );
            }
        }

        StructuredGrid gs_global( gridGlobal_ );
        ATLAS_ASSERT( gs_global );  // assert structured grid
        StructuredGrid gsLeg = ( useGlobalLeg ? gs_global : g );
        nlonsMaxGlobal_      = gs_global.nxmax();
        jlonMin_.resize( 1 );
        jlonMin_[0]  = 0;
        jlatMin_     = 0;
        nlatsGlobal_ = gs_global.ny();
        if ( grid_.domain().global() ) { Log::debug() << "Global grid with " << nlats << " latitudes." << std::endl; }
        else {
            Log::debug() << "Grid has " << nlats << " latitudes. Global grid has " << nlatsGlobal_ << std::endl;
        }
        if ( useGlobalLeg ) { nlatsLeg_ = ( nlatsGlobal_ + 1 ) / 2; }
        else {
            nlatsLeg_        = nlatsLegDomain_;
            nlatsLegReduced_ = nlatsLeg_;
        }
        for ( int jlat = 0; jlat < nlatsGlobal_; jlat++ ) {
            if ( gs_global.y( jlat ) > g.y( 0 ) ) {
                //Log::info() << gs_global.y( jlat ) << ">" << g.y( 0 ) << " ";
                jlatMin_++;
            };
        }
        //Log::info() << std::endl;
        int jlatMinLeg_ = jlatMin_;
        if ( nlatsNH_ < nlatsSH_ ) { jlatMinLeg_ += nlatsNH_ - nlatsSH_; };
        if ( jlatMin_ >= ( nlatsGlobal_ + 1 ) / 2 ) {
            jlatMinLeg_ -= 2 * ( jlatMin_ - ( nlatsGlobal_ + 1 ) / 2 );
            if ( nlatsGlobal_ % 2 == 1 ) { jlatMinLeg_--; }
        };
        if ( useGlobalLeg ) { nlatsLegReduced_ = jlatMinLeg_ + nlatsLegDomain_; }

        // reduce truncation towards the pole for reduced meshes:
        nlat0_.resize( truncation_ + 1 );
        if ( no_nest ) {
            for ( int j = 0; j <= truncation_; j++ ) {
                nlat0_[j] = 0;
            }
        }
        else {
            int nmen0 = -1;
            for ( int jlat = 0; jlat < nlatsGlobal_ / 2; jlat++ ) {
                double lat = gs_global.y( jlat ) * util::Constants::degreesToRadians();
                int nmen = fourier_truncation( truncation_, gs_global.nx( jlat ), gs_global.nxmax(), nlatsGlobal_, lat,
                                               RegularGrid( gs_global ) );
                nmen     = std::max( nmen0, nmen );
                int ndgluj = nlatsLeg_ - std::min( nlatsLeg_, nlatsLeg_ + jlatMinLeg_ - jlat );
                if ( useGlobalLeg ) { ndgluj = std::max( jlatMinLeg_, jlat ); }
                for ( int j = nmen0 + 1; j <= nmen; j++ ) {
                    nlat0_[j] = ndgluj;
                }
                nmen0 = nmen;
            }
            for ( int j = nmen0 + 1; j <= truncation_; j++ ) {
                nlat0_[j] = nlatsLeg_;
            }
        }
        /*Log::info() << "nlats=" << g.ny() << " nlatsGlobal=" << gs_global.ny() << " jlatMin=" << jlatMin_
                    << " jlatMinLeg=" << jlatMinLeg_ << " nlatsGlobal/2-nlatsLeg=" << nlatsGlobal_ / 2 - nlatsLeg_
                    << " nlatsLeg_=" << nlatsLeg_ << " nlatsLegDomain_=" << nlatsLegDomain_ << std::endl;*/

        // compute longitudinal location of domain within global grid for using FFT:
        auto wrapAngle = [&]( double angle ) {
            double result = std::fmod( angle, 360. );
            if ( result < 0. ) { result += 360.; }
            return result;
        };
        if ( useFFT_ ) {
            double lonmin = wrapAngle( g.x( 0, 0 ) );
            if ( nlonsMax < fft_threshold * nlonsMaxGlobal_ ) { useFFT_ = false; }
            else {
                // need to use FFT with cropped grid
                if ( RegularGrid( gridGlobal_ ) ) {
                    for ( idx_t jlon = 0; jlon < nlonsMaxGlobal_; ++jlon ) {
                        if ( gs_global.x( jlon, 0 ) < lonmin ) { jlonMin_[0]++; }
                    }
                }
                else {
                    nlonsGlobal_.resize( nlats );
                    jlonMin_.resize( nlats );
                    for ( idx_t jlat = 0; jlat < nlats; jlat++ ) {
                        double lonmin      = wrapAngle( g.x( 0, jlat ) );
                        nlonsGlobal_[jlat] = gs_global.nx( jlat + jlatMin_ );
                        jlonMin_[jlat]     = 0;
                        for ( idx_t jlon = 0; jlon < nlonsGlobal_[jlat]; ++jlon ) {
                            if ( gs_global.x( jlon, jlat + jlatMin_ ) < lonmin ) { jlonMin_[jlat]++; }
                        }
                    }
                }
            }
        }
        //Log::info() << "nlats=" << g.ny() << " nlatsGlobal=" << gs_global.ny() << std::endl;
        std::vector<double> lats( nlatsLeg_ );
        std::vector<double> lons( nlonsMax );
        if ( nlatsNH_ >= nlatsSH_ || useGlobalLeg ) {
            for ( idx_t j = 0; j < nlatsLeg_; ++j ) {
                double lat = gsLeg.y( j );
                if ( lat > latPole ) { lat = latPole; }
                if ( lat < -latPole ) { lat = -latPole; }
                lats[j] = lat * util::Constants::degreesToRadians();
            }
        }
        else {
            for ( idx_t j = nlats - 1, idx = 0; idx < nlatsLeg_; --j, ++idx ) {
                double lat = gsLeg.y( j );
                if ( lat > latPole ) { lat = latPole; }
                if ( lat < -latPole ) { lat = -latPole; }
                lats[idx] = -lat * util::Constants::degreesToRadians();
            }
        }
        for ( idx_t j = 0; j < nlonsMax; ++j ) {
            lons[j] = g.x( j, 0 ) * util::Constants::degreesToRadians();
        }
        /*Log::info() << "lats: ";
        for ( int j = 0; j < nlatsLeg_; j++ ) {
            Log::info() << lats[j] << " ";
        }
        Log::info() << std::endl;*/

        // precomputations for Legendre polynomials:
        {
            const auto nlatsLeg = size_t( nlatsLeg_ );
            size_t size_sym     = 0;
            size_t size_asym    = 0;
            legendre_sym_begin_.resize( truncation_ + 3 );
            legendre_asym_begin_.resize( truncation_ + 3 );
            legendre_sym_begin_[0]  = 0;
            legendre_asym_begin_[0] = 0;
            for ( idx_t jm = 0; jm <= truncation_ + 1; jm++ ) {
                size_sym += add_padding( num_n( truncation_ + 1, jm, /*symmetric*/ true ) * nlatsLeg );
                size_asym += add_padding( num_n( truncation_ + 1, jm, /*symmetric*/ false ) * nlatsLeg );
                legendre_sym_begin_[jm + 1]  = size_sym;
                legendre_asym_begin_[jm + 1] = size_asym;
            }

            if ( legendre_cache_ ) {
                ReadCache legendre( legendre_cache_ );
                legendre_sym_  = legendre.read<double>( size_sym );
                legendre_asym_ = legendre.read<double>( size_asym );
                ATLAS_ASSERT( legendre.pos == legendre_cachesize_ );
                // TODO: check this is all aligned...
            }
            else {
                if ( TransParameters( config ).export_legendre() ) {
                    ATLAS_ASSERT( not cache_.legendre() );

                    size_t bytes = sizeof( double ) * ( size_sym + size_asym );
                    Log::debug() << "TransLocal: allocating LegendreCache: " << eckit::Bytes( bytes ) << std::endl;
                    export_legendre_ = LegendreCache( bytes );

                    legendre_cachesize_ = export_legendre_.legendre().size();
                    legendre_cache_     = export_legendre_.legendre().data();
                    ReadCache legendre( legendre_cache_ );
                    legendre_sym_  = legendre.read<double>( size_sym );
                    legendre_asym_ = legendre.read<double>( size_asym );
                }
                else {
                    alloc_aligned( legendre_sym_, size_sym, "symmetric" );
                    alloc_aligned( legendre_asym_, size_asym, "asymmetric" );
                }

                ATLAS_TRACE_SCOPE( "Legendre precomputations (structured)" ) {
                    compute_legendre_polynomials( truncation_ + 1, nlatsLeg_, lats.data(), legendre_sym_,
                                                  legendre_asym_, legendre_sym_begin_.data(),
                                                  legendre_asym_begin_.data() );
                }
                std::string file_path = TransParameters( config ).write_legendre();
                if ( file_path.size() ) {
                    ATLAS_TRACE( "Write LegendreCache to file" );
                    Log::debug() << "Writing Legendre cache file ..." << std::endl;
                    Log::debug() << "    path: " << file_path << std::endl;
                    WriteCache legendre( file_path );
                    legendre.write( legendre_sym_, size_sym );
                    legendre.write( legendre_asym_, size_asym );
                    Log::debug() << "    size: " << eckit::Bytes( legendre.pos ) << std::endl;
                }
            }
        }

        // precomputations for Fourier transformations:
        if ( useFFT_ ) {
#if ATLAS_HAVE_FFTW && !TRANSLOCAL_DGEMM2
            {
                ATLAS_TRACE( "Fourier precomputations (FFTW)" );
                int num_complex = ( nlonsMaxGlobal_ / 2 ) + 1;
                fftw_->in       = fftw_alloc_complex( nlats * num_complex );
                fftw_->out      = fftw_alloc_real( nlats * nlonsMaxGlobal_ );

                if ( fft_cache_ ) {
                    Log::debug() << "Import FFTW wisdom from cache" << std::endl;
                    fftw_import_wisdom_from_string( static_cast<const char*>( fft_cache_ ) );
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
                if ( RegularGrid( gridGlobal_ ) ) {
                    fftw_->plans.resize( 1 );
                    fftw_->plans[0] =
                        fftw_plan_many_dft_c2r( 1, &nlonsMaxGlobal_, nlats, fftw_->in, nullptr, 1, num_complex,
                                                fftw_->out, nullptr, 1, nlonsMaxGlobal_, FFTW_ESTIMATE );
                }
                else {
                    fftw_->plans.resize( nlatsLegDomain_ );
                    for ( int j = 0; j < nlatsLegDomain_; j++ ) {
                        int nlonsGlobalj = gs_global.nx( jlatMinLeg_ + j );
                        //ASSERT( nlonsGlobalj > 0 && nlonsGlobalj <= nlonsMaxGlobal_ );
                        fftw_->plans[j] = fftw_plan_dft_c2r_1d( nlonsGlobalj, fftw_->in, fftw_->out, FFTW_ESTIMATE );
                    }
                }
                std::string file_path = TransParameters( config ).write_fft();
                if ( file_path.size() ) {
                    Log::debug() << "Write FFTW wisdom to file " << file_path << std::endl;
                    //bool success = fftw_export_wisdom_to_filename( "wisdom.bin" );
                    //ASSERT( success );
                    //std::ofstream write( file_path );
                    //write << FFTW_Wisdom();

                    FILE* file_fftw = fopen( file_path.c_str(), "wb" );
                    fftw_export_wisdom_to_file( file_fftw );
                    fclose( file_fftw );
                }
                //                std::string newWisdom( fftw_export_wisdom_to_string() );
                //                if ( 1.1 * wisdomString.length() < newWisdom.length() ) {
                //                    std::ofstream write( "wisdom.bin" );
                //                    write << newWisdom;
                //                    write.close();
                //                }
            }
                // other FFT implementations should be added with #elif statements
#else
            useFFT_               = false;  // no FFT implemented => default to dgemm
            std::string file_path = TransParameters( config ).write_fft();
            if ( file_path.size() ) {
                std::ofstream write( file_path );
                write << "No cache available, as FFTW is not enabled" << std::endl;
                write.close();
            }

#endif
        }
        if ( !useFFT_ ) {
            Log::warning()
                << "WARNING: Spectral transform results may contain aliasing errors. This will be addressed soon."
                << std::endl;

            alloc_aligned( fourier_, 2 * ( truncation_ + 1 ) * nlonsMax, "Fourier coeffs." );
#if !TRANSLOCAL_DGEMM2
            {
                ATLAS_TRACE( "Fourier precomputations (NoFFT)" );
                int idx = 0;
                for ( int jm = 0; jm < truncation_ + 1; jm++ ) {
                    double factor = 1.;
                    if ( jm > 0 ) { factor = 2.; }
                    for ( int jlon = 0; jlon < nlonsMax; jlon++ ) {
                        fourier_[idx++] = +std::cos( jm * lons[jlon] ) * factor;  // real part
                    }
                    for ( int jlon = 0; jlon < nlonsMax; jlon++ ) {
                        fourier_[idx++] = -std::sin( jm * lons[jlon] ) * factor;  // imaginary part
                    }
                }
            }
#else
            {
                ATLAS_TRACE( "precomp Fourier" );
                int idx = 0;
                for ( int jlon = 0; jlon < nlonsMax; jlon++ ) {
                    double factor = 1.;
                    for ( int jm = 0; jm < truncation_ + 1; jm++ ) {
                        if ( jm > 0 ) { factor = 2.; }
                        fourier_[idx++] = +std::cos( jm * lons[jlon] ) * factor;  // real part
                        fourier_[idx++] = -std::sin( jm * lons[jlon] ) * factor;  // imaginary part
                    }
                }
            }
#endif
        }
    }
    else {
        // unstructured grid
        if ( unstruct_precomp_ ) {
            ATLAS_TRACE( "Legendre precomputations (unstructured)" );

            if ( warning() ) {
                Log::warning()
                    << "WARNING: Precomputations for spectral transforms could take a long time as there's no structure"
                       " to take advantage of!!!"
                    << std::endl
                    << "The precomputed values consume at least "
                    << eckit::Bytes( sizeof( double ) * legendre_size( truncation_ ) * grid_.size() ) << " ("
                    << eckit::Bytes( sizeof( double ) * legendre_size( truncation_ ) ) << " for each of "
                    << grid_.size() << " grid points )" << std::endl
                    << "Furthermore, results may contain aliasing errors." << std::endl;
            }

            std::vector<double> lats( grid_.size() );
            alloc_aligned( legendre_, legendre_size( truncation_ ) * grid_.size(), "Legendre coeffs." );
            int j( 0 );
            for ( PointLonLat p : grid_.lonlat() ) {
                lats[j++] = p.lat() * util::Constants::degreesToRadians();
            }
            compute_legendre_polynomials_all( truncation_, grid_.size(), lats.data(), legendre_ );
        }
        if ( TransParameters( config ).write_legendre().size() ) {
            throw_NotImplemented(
                "Caching for unstructured grids or structured grids with projections not yet implemented", Here() );
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

TransLocal::TransLocal( const Grid& grid, const long truncation, const eckit::Configuration& config ) :
    TransLocal( Cache(), grid, grid.domain(), truncation, config ) {}

TransLocal::TransLocal( const Grid& grid, const Domain& domain, const long truncation,
                        const eckit::Configuration& config ) :
    TransLocal( Cache(), grid, domain, truncation, config ) {}

TransLocal::TransLocal( const Cache& cache, const Grid& grid, const long truncation,
                        const eckit::Configuration& config ) :
    TransLocal( cache, grid, grid.domain(), truncation, config ) {}

// --------------------------------------------------------------------------------------------------------------------

TransLocal::~TransLocal() {
    if ( StructuredGrid( grid_ ) && not grid_.projection() ) {
        if ( not legendre_cache_ ) {
            free_aligned( legendre_sym_, "symmetric" );
            free_aligned( legendre_asym_, "asymmetric" );
        }
        if ( useFFT_ ) {
#if ATLAS_HAVE_FFTW && !TRANSLOCAL_DGEMM2
            for ( idx_t j = 0, size = static_cast<idx_t>( fftw_->plans.size() ); j < size; j++ ) {
                fftw_destroy_plan( fftw_->plans[j] );
            }
            fftw_free( fftw_->in );
            fftw_free( fftw_->out );
#endif
        }
        else {
            free_aligned( fourier_, "Fourier coeffs." );
        }
    }
    else {
        if ( unstruct_precomp_ ) { free_aligned( legendre_, "Legendre coeffs." ); }
    }
}

// --------------------------------------------------------------------------------------------------------------------

const functionspace::Spectral& TransLocal::spectral() const {
    if ( not spectral_ ) { spectral_ = functionspace::Spectral( Trans( this ) ); }
    return spectral_;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans( const Field& spfield, Field& gpfield, const eckit::Configuration& config ) const {

// VERY PRELIMINARY IMPLEMENTATION WITHOUT ANY GUARANTEES    
    int nb_scalar_fields = 1;
    const auto scalar_spectra = array::make_view<double,1>( spfield );
    auto gp_fields = array::make_view<double,1>( spfield );

    if( gp_fields.shape(0) != grid().size() ) {
        ATLAS_DEBUG_VAR( gp_fields.shape(0) );
        ATLAS_DEBUG_VAR( grid().size() );
        ATLAS_ASSERT( gp_fields.shape(0) == grid().size() );
    }

    invtrans( nb_scalar_fields, scalar_spectra.data(), gp_fields.data(), config );

}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans( const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& config ) const {

// VERY PRELIMINARY IMPLEMENTATION WITHOUT ANY GUARANTEES    
    ATLAS_ASSERT( spfields.size() == gpfields.size() );
    for( idx_t f=0; f<spfields.size(); ++f ) {
        invtrans( spfields[f], gpfields[f], config );
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_grad( const Field& /*spfield*/, Field& /*gradfield*/, const eckit::Configuration& ) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_grad( const FieldSet& /*spfields*/, FieldSet& /*gradfields*/,
                                const eckit::Configuration& ) const {
    ATLAS_NOTIMPLEMENTED;
}

// --------------------------------------------------------------------------------------------------------------------

void gp_transpose( const int nb_size, const int nb_fields, const double gp_tmp[], double gp_fields[] ) {
    for ( int jgp = 0; jgp < nb_size; jgp++ ) {
        for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
            gp_fields[jfld * nb_size + jgp] = gp_tmp[jgp * nb_fields + jfld];
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_vordiv2wind( const Field& spvor, const Field& spdiv, Field& gpwind,
                                       const eckit::Configuration& config ) const {

// VERY PRELIMINARY IMPLEMENTATION WITHOUT ANY GUARANTEES    
    int nb_vordiv_fields = 1;
    const auto vorticity_spectra = array::make_view<double,1>( spvor );
    const auto divergence_spectra = array::make_view<double,1>( spdiv );
    auto gp_fields = array::make_view<double,2>( gpwind );

    if( gp_fields.shape(1) == grid().size() && gp_fields.shape(0) == 2 ) {
        invtrans( nb_vordiv_fields, vorticity_spectra.data(), divergence_spectra.data(), gp_fields.data(), config );
    }
    else if( gp_fields.shape(0) == grid().size() && gp_fields.shape(1) == 2 ) {
        array::ArrayT<double> gpwind_t( gp_fields.shape(1), gp_fields.shape(0) );
        auto gp_fields_t = array::make_view<double,2>( gpwind_t );
        invtrans( nb_vordiv_fields, vorticity_spectra.data(), divergence_spectra.data(), gp_fields_t.data(), config );
        gp_transpose( grid().size(), 2, gp_fields_t.data(), gp_fields.data() );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans( const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                           const eckit::Configuration& config ) const {
    invtrans_uv( truncation_, nb_scalar_fields, 0, scalar_spectra, gp_fields, config );
}


// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_legendre( const int truncation, const int nlats, const int nb_fields,
                                    const int /*nb_vordiv_fields*/, const double scalar_spectra[], double scl_fourier[],
                                    const eckit::Configuration& ) const {
    // Legendre transform:
    {
        Log::debug() << "Legendre dgemm: using " << nlatsLegReduced_ - nlat0_[0] << " latitudes out of "
                     << nlatsGlobal_ / 2 << std::endl;
        ATLAS_TRACE( "Inverse Legendre Transform (GEMM)" );
        for ( int jm = 0; jm <= truncation_; jm++ ) {
            size_t size_sym  = num_n( truncation_ + 1, jm, true );
            size_t size_asym = num_n( truncation_ + 1, jm, false );
            const int n_imag = ( jm ? 2 : 1 );
            int size_fourier = nb_fields * n_imag * ( nlatsLegReduced_ - nlat0_[jm] );
            if ( size_fourier > 0 ) {
                auto posFourier = [&]( int jfld, int imag, int jlat, int jm, int nlatsH ) {
                    return jfld + nb_fields * ( imag + n_imag * ( nlatsLegReduced_ - nlat0_[jm] - nlatsH + jlat ) );
                };
                double* scalar_sym;
                double* scalar_asym;
                double* scl_fourier_sym;
                double* scl_fourier_asym;
                alloc_aligned( scalar_sym, n_imag * nb_fields * size_sym );
                alloc_aligned( scalar_asym, n_imag * nb_fields * size_asym );
                alloc_aligned( scl_fourier_sym, size_fourier );
                alloc_aligned( scl_fourier_asym, size_fourier );
                {
                    //ATLAS_TRACE( "Legendre split" );
                    idx_t idx = 0, is = 0, ia = 0, ioff = ( 2 * truncation + 3 - jm ) * jm / 2 * nb_fields * 2;
                    // the choice between the following two code lines determines whether
                    // total wavenumbers are summed in an ascending or descending order.
                    // The trans library in IFS uses descending order because it should
                    // be more accurate (higher wavenumbers have smaller contributions).
                    // This also needs to be changed when splitting the spectral data in
                    // compute_legendre_polynomials!
                    //for ( int jn = jm; jn <= truncation_ + 1; jn++ ) {
                    for ( int jn = truncation_ + 1; jn >= jm; jn-- ) {
                        for ( int imag = 0; imag < n_imag; imag++ ) {
                            for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                                idx = jfld + nb_fields * ( imag + 2 * ( jn - jm ) );
                                if ( jn <= truncation && jm < truncation ) {
                                    if ( ( jn - jm ) % 2 == 0 ) { scalar_sym[is++] = scalar_spectra[idx + ioff]; }
                                    else {
                                        scalar_asym[ia++] = scalar_spectra[idx + ioff];
                                    }
                                }
                                else {
                                    if ( ( jn - jm ) % 2 == 0 ) { scalar_sym[is++] = 0.; }
                                    else {
                                        scalar_asym[ia++] = 0.;
                                    }
                                }
                            }
                        }
                    }
                    ATLAS_ASSERT( size_t( ia ) == n_imag * nb_fields * size_asym &&
                                  size_t( is ) == n_imag * nb_fields * size_sym );
                }
                if ( nlatsLegReduced_ - nlat0_[jm] > 0 ) {
                    {
                        eckit::linalg::Matrix A( scalar_sym, nb_fields * n_imag, size_sym );
                        eckit::linalg::Matrix B( legendre_sym_ + legendre_sym_begin_[jm] + nlat0_[jm] * size_sym,
                                                 size_sym, nlatsLegReduced_ - nlat0_[jm] );
                        eckit::linalg::Matrix C( scl_fourier_sym, nb_fields * n_imag, nlatsLegReduced_ - nlat0_[jm] );
                        linalg_.gemm( A, B, C );
                        /*Log::info() << "sym: ";
                        for ( int j = 0; j < size_sym * ( nlatsLegReduced_ - nlat0_[jm] ); j++ ) {
                            Log::info() << legendre_sym_[j + legendre_sym_begin_[jm] + nlat0_[jm] * size_sym] << " ";
                        }
                        Log::info() << std::endl;*/
                    }
                    if ( size_asym > 0 ) {
                        eckit::linalg::Matrix A( scalar_asym, nb_fields * n_imag, size_asym );
                        eckit::linalg::Matrix B( legendre_asym_ + legendre_asym_begin_[jm] + nlat0_[jm] * size_asym,
                                                 size_asym, nlatsLegReduced_ - nlat0_[jm] );
                        eckit::linalg::Matrix C( scl_fourier_asym, nb_fields * n_imag, nlatsLegReduced_ - nlat0_[jm] );
                        linalg_.gemm( A, B, C );
                        /*Log::info() << "asym: ";
                        for ( int j = 0; j < size_asym * ( nlatsLegReduced_ - nlat0_[jm] ); j++ ) {
                            Log::info() << legendre_asym_[j + legendre_asym_begin_[jm] + nlat0_[jm] * size_asym] << " ";
                        }
                        Log::info() << std::endl;*/
                    }
                }
                {
                    //ATLAS_TRACE( "merge spheres" );
                    // northern hemisphere:
                    for ( int jlat = 0; jlat < nlatsNH_; jlat++ ) {
                        if ( nlatsLegReduced_ - nlat0_[jm] - nlatsNH_ + jlat >= 0 ) {
                            for ( int imag = 0; imag < n_imag; imag++ ) {
                                for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                                    int idx = posFourier( jfld, imag, jlat, jm, nlatsNH_ );
                                    scl_fourier[posMethod( jfld, imag, jlat, jm, nb_fields, nlats )] =
                                        scl_fourier_sym[idx] + scl_fourier_asym[idx];
                                }
                            }
                        }
                        else {
                            for ( int imag = 0; imag < n_imag; imag++ ) {
                                for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                                    scl_fourier[posMethod( jfld, imag, jlat, jm, nb_fields, nlats )] = 0.;
                                }
                            }
                        }
                        /*for ( int imag = 0; imag < n_imag; imag++ ) {
                        for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                            if ( scl_fourier[posMethod( jfld, imag, jlat, jm, nb_fields, nlats )] > 0. ) {
                                Log::info() << "jm=" << jm << " jlat=" << jlat << " nlatsLeg_=" << nlatsLeg_
                                            << " nlat0=" << nlat0_[jm] << " nlatsNH=" << nlatsNH_ << std::endl;
                            }
                        }
                    }*/
                    }
                    // southern hemisphere:
                    for ( int jlat = 0; jlat < nlatsSH_; jlat++ ) {
                        int jslat = nlats - jlat - 1;
                        if ( nlatsLegReduced_ - nlat0_[jm] - nlatsSH_ + jlat >= 0 ) {
                            for ( int imag = 0; imag < n_imag; imag++ ) {
                                for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                                    int idx = posFourier( jfld, imag, jlat, jm, nlatsSH_ );
                                    scl_fourier[posMethod( jfld, imag, jslat, jm, nb_fields, nlats )] =
                                        scl_fourier_sym[idx] - scl_fourier_asym[idx];
                                }
                            }
                        }
                        else {
                            for ( int imag = 0; imag < n_imag; imag++ ) {
                                for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                                    scl_fourier[posMethod( jfld, imag, jslat, jm, nb_fields, nlats )] = 0.;
                                }
                            }
                        }
                    }
                }
                free_aligned( scalar_sym );
                free_aligned( scalar_asym );
                free_aligned( scl_fourier_sym );
                free_aligned( scl_fourier_asym );
            }
            else {
                for ( int jlat = 0; jlat < nlats; jlat++ ) {
                    for ( int imag = 0; imag < n_imag; imag++ ) {
                        for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                            scl_fourier[posMethod( jfld, imag, jlat, jm, nb_fields, nlats )] = 0.;
                        }
                    }
                }
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_fourier_regular( const int nlats, const int nlons, const int nb_fields, double scl_fourier[],
                                           double gp_fields[], const eckit::Configuration& ) const {
    // Fourier transformation:
    if ( useFFT_ ) {
#if ATLAS_HAVE_FFTW && !TRANSLOCAL_DGEMM2
        {
            int num_complex = ( nlonsMaxGlobal_ / 2 ) + 1;
            {
                ATLAS_TRACE( "Inverse Fourier Transform (FFTW, RegularGrid)" );
                for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                    int idx = 0;
                    for ( int jlat = 0; jlat < nlats; jlat++ ) {
                        fftw_->in[idx++][0] = scl_fourier[posMethod( jfld, 0, jlat, 0, nb_fields, nlats )];
                        for ( int jm = 1; jm < num_complex; jm++, idx++ ) {
                            for ( int imag = 0; imag < 2; imag++ ) {
                                if ( jm <= truncation_ ) {
                                    fftw_->in[idx][imag] =
                                        scl_fourier[posMethod( jfld, imag, jlat, jm, nb_fields, nlats )];
                                }
                                else {
                                    fftw_->in[idx][imag] = 0.;
                                }
                            }
                        }
                    }
                    fftw_execute_dft_c2r( fftw_->plans[0], fftw_->in, fftw_->out );
                    for ( int jlat = 0; jlat < nlats; jlat++ ) {
                        for ( int jlon = 0; jlon < nlons; jlon++ ) {
                            int j = jlon + jlonMin_[0];
                            if ( j >= nlonsMaxGlobal_ ) { j -= nlonsMaxGlobal_; }
                            gp_fields[jlon + nlons * ( jlat + nlats * jfld )] = fftw_->out[j + nlonsMaxGlobal_ * jlat];
                        }
                    }
                }
            }
        }
#endif
    }
    else {
#if !TRANSLOCAL_DGEMM2
        // dgemm-method 1
        {
            ATLAS_TRACE( "Inverse Fourier Transform (NoFFT)" );
            eckit::linalg::Matrix A( fourier_, nlons, ( truncation_ + 1 ) * 2 );
            eckit::linalg::Matrix B( scl_fourier, ( truncation_ + 1 ) * 2, nb_fields * nlats );
            eckit::linalg::Matrix C( gp_fields, nlons, nb_fields * nlats );

            linalg_.gemm( A, B, C );
        }
#else
        // dgemm-method 2
        // should be faster for small domains or large truncation
        // but have not found any significant speedup so far
        double* gp;
        alloc_aligned( gp, nb_fields * grid_.size() );
        {
            ATLAS_TRACE( "Fourier dgemm method 2" );
            eckit::linalg::Matrix A( scl_fourier, nb_fields * nlats, ( truncation_ + 1 ) * 2 );
            eckit::linalg::Matrix B( fourier_, ( truncation_ + 1 ) * 2, nlons );
            eckit::linalg::Matrix C( gp, nb_fields * nlats, nlons );
            linalg_.gemm( A, B, C );
        }

        // Transposition in grid point space:
        {
            ATLAS_TRACE( "transposition in gp-space" );
            int idx = 0;
            for ( int jlon = 0; jlon < nlons; jlon++ ) {
                for ( int jlat = 0; jlat < nlats; jlat++ ) {
                    for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                        int pos_tp = jlon + nlons * ( jlat + nlats * ( jfld ) );
                        //int pos  = jfld + nb_fields * ( jlat + nlats * ( jlon ) );
                        gp_fields[pos_tp] = gp[idx++];  // = gp[pos]
                    }
                }
            }
        }
        free_aligned( gp );
#endif
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_fourier_reduced( const int nlats, const StructuredGrid& g, const int nb_fields,
                                           double scl_fourier[], double gp_fields[],
                                           const eckit::Configuration& ) const {
    // Fourier transformation:
    if ( useFFT_ ) {
#if ATLAS_HAVE_FFTW && !TRANSLOCAL_DGEMM2
        {
            {
                ATLAS_TRACE( "Inverse Fourier Transform (FFTW, ReducedGrid)" );
                int jgp = 0;
                for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                    for ( int jlat = 0; jlat < nlats; jlat++ ) {
                        int idx = 0;
                        //Log::info() << jlat << "in:" << std::endl;
                        int num_complex     = ( nlonsGlobal_[jlat] / 2 ) + 1;
                        fftw_->in[idx++][0] = scl_fourier[posMethod( jfld, 0, jlat, 0, nb_fields, nlats )];
                        //Log::info() << fftw_->in[0][0] << " ";
                        for ( int jm = 1; jm < num_complex; jm++, idx++ ) {
                            for ( int imag = 0; imag < 2; imag++ ) {
                                if ( jm <= truncation_ ) {
                                    fftw_->in[idx][imag] =
                                        scl_fourier[posMethod( jfld, imag, jlat, jm, nb_fields, nlats )];
                                }
                                else {
                                    fftw_->in[idx][imag] = 0.;
                                }
                                //Log::info() << fftw_->in[idx][imag] << " ";
                            }
                        }
                        //Log::info() << std::endl;
                        //Log::info() << jlat << "out:" << std::endl;
                        int jplan = nlatsLegDomain_ - nlatsNH_ + jlat;
                        if ( jplan >= nlatsLegDomain_ ) { jplan = nlats - 1 + nlatsLegDomain_ - nlatsSH_ - jlat; };
                        //ASSERT( jplan < nlatsLeg_ && jplan >= 0 );
                        fftw_execute_dft_c2r( fftw_->plans[jplan], fftw_->in, fftw_->out );
                        for ( int jlon = 0; jlon < g.nx( jlat ); jlon++ ) {
                            int j = jlon + jlonMin_[jlat];
                            if ( j >= nlonsGlobal_[jlat] ) { j -= nlonsGlobal_[jlat]; }
                            //Log::info() << fftw_->out[j] << " ";
                            ATLAS_ASSERT( j < nlonsMaxGlobal_ );
                            gp_fields[jgp++] = fftw_->out[j];
                        }
                        //Log::info() << std::endl;
                    }
                }
            }
        }
#endif
    }
    else {
        throw_NotImplemented(
            "Using dgemm in Fourier transform for reduced grids is extremely slow. Please install and use FFTW!",
            Here() );
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_unstructured_precomp( const int truncation, const int nb_fields, const int nb_vordiv_fields,
                                                const double scalar_spectra[], double gp_fields[],
                                                const eckit::Configuration& ) const {
    ATLAS_TRACE( "invtrans_uv unstructured" );

    const int nlats        = grid_.size();
    const int size_fourier = nb_fields * 2;
    double* scl_fourier;
    double* scl_fourier_tp;
    double* fouriertp;
    double* gp_opt;
    alloc_aligned( scl_fourier, size_fourier * (truncation)*nlats );
    alloc_aligned( scl_fourier_tp, size_fourier * ( truncation ) );
    alloc_aligned( fouriertp, 2 * ( truncation ) );
    alloc_aligned( gp_opt, nb_fields );

    {
        ATLAS_TRACE( "Inverse Legendre Transform (GEMM)" );
        for ( int jm = 0; jm < truncation; jm++ ) {
            const int noff = ( 2 * truncation + 3 - jm ) * jm / 2, ns = truncation - jm + 1;
            eckit::linalg::Matrix A( eckit::linalg::Matrix(
                const_cast<double*>( scalar_spectra ) + nb_fields * 2 * noff, nb_fields * 2, ns ) );
            eckit::linalg::Matrix B( legendre_ + noff * nlats, ns, nlats );
            eckit::linalg::Matrix C( scl_fourier + jm * size_fourier * nlats, nb_fields * 2, nlats );
            linalg_.gemm( A, B, C );
        }
    }

    // loop over all points:
    {
        ATLAS_TRACE( "Inverse Fourier Transform (NoFFT)" );
        int ip = 0;
        for ( const PointLonLat p : grid_.lonlat() ) {
            const double lon = p.lon() * util::Constants::degreesToRadians();
            const double lat = p.lat() * util::Constants::degreesToRadians();
            {
                //ATLAS_TRACE( "opt transposition in Fourier" );
                for ( int jm = 0; jm < truncation; jm++ ) {
                    int idx = nb_fields * 2 * ( ip + nlats * jm );
                    for ( int imag = 0; imag < 2; imag++ ) {
                        for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                            const int pos_tp = imag + 2 * ( jm + ( truncation ) * ( jfld ) );
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
                for ( int jm = 1; jm < truncation; jm++ ) {
                    fouriertp[idx++] = +2. * std::cos( jm * lon );  // real part
                    fouriertp[idx++] = -2. * std::sin( jm * lon );  // imaginary part
                }
            }
            {
                //ATLAS_TRACE( "opt Fourier dgemm" );
                eckit::linalg::Matrix A( fouriertp, 1, (truncation)*2 );
                eckit::linalg::Matrix B( scl_fourier_tp, (truncation)*2, nb_fields );
                eckit::linalg::Matrix C( gp_opt, 1, nb_fields );
                linalg_.gemm( A, B, C );
                for ( int j = 0; j < nb_fields; j++ ) {
                    gp_fields[ip + j * grid_.size()] = gp_opt[j];
                }
            }
            // Computing u,v from U,V:
            {
                if ( nb_vordiv_fields > 0 ) {
                    //ATLAS_TRACE( " u,v from U,V" );
                    double coslat = std::cos( lat );
                    for ( int j = 0; j < 2 * nb_vordiv_fields && j < nb_fields; j++ ) {
                        gp_fields[ip + j * grid_.size()] /= coslat;
                    }
                }
            }
            ++ip;
        }
    }
    free_aligned( scl_fourier );
    free_aligned( scl_fourier_tp );
    free_aligned( fouriertp );
    free_aligned( gp_opt );
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_unstructured( const int truncation, const int nb_fields, const int nb_vordiv_fields,
                                        const double scalar_spectra[], double gp_fields[],
                                        const eckit::Configuration& config ) const {
    ATLAS_TRACE( "invtrans_unstructured" );

    if ( warning( config ) ) {
        Log::warning() << "WARNING: Spectral transforms could take a long time (unstructured grid approach). Results "
                          "may contain aliasing errors."
                       << std::endl;
    }

    double* zfn;
    alloc_aligned( zfn, ( truncation + 1 ) * ( truncation + 1 ) );
    compute_zfn( truncation, zfn );
    int size_fourier = nb_fields * 2;
    double* legendre;
    double* scl_fourier;
    double* scl_fourier_tp;
    double* fouriertp;
    double* gp_opt;
    alloc_aligned( legendre, legendre_size( truncation + 1 ) );
    alloc_aligned( scl_fourier, size_fourier * ( truncation + 1 ) );
    alloc_aligned( scl_fourier_tp, size_fourier * ( truncation + 1 ) );
    alloc_aligned( fouriertp, 2 * ( truncation + 1 ) );
    alloc_aligned( gp_opt, nb_fields );


    // loop over all points:
    int ip = 0;
    for ( const PointLonLat p : grid_.lonlat() ) {
        const double lon = p.lon() * util::Constants::degreesToRadians();
        const double lat = p.lat() * util::Constants::degreesToRadians();
        compute_legendre_polynomials_lat( truncation, lat, legendre, zfn );
        // Legendre transform:
        {
            //ATLAS_TRACE( "opt Legendre dgemm" );
            for ( int jm = 0; jm <= truncation; jm++ ) {
                const int noff = ( 2 * truncation + 3 - jm ) * jm / 2, ns = truncation - jm + 1;
                eckit::linalg::Matrix A( eckit::linalg::Matrix(
                    const_cast<double*>( scalar_spectra ) + nb_fields * 2 * noff, nb_fields * 2, ns ) );
                eckit::linalg::Matrix B( legendre + noff, ns, 1 );
                eckit::linalg::Matrix C( scl_fourier + jm * size_fourier, nb_fields * 2, 1 );
                linalg_.gemm( A, B, C );
            }
        }
        {
            //ATLAS_TRACE( "opt transposition in Fourier" );
            int idx = 0;
            for ( int jm = 0; jm < truncation + 1; jm++ ) {
                for ( int imag = 0; imag < 2; imag++ ) {
                    for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                        const int pos_tp = imag + 2 * ( jm + ( truncation + 1 ) * ( jfld ) );
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
        for ( int jm = 1; jm < truncation + 1; jm++ ) {
            fouriertp[idx++] = +2. * std::cos( jm * lon );  // real part
            fouriertp[idx++] = -2. * std::sin( jm * lon );  // imaginary part
        }
        {
            //ATLAS_TRACE( "opt Fourier dgemm" );
            eckit::linalg::Matrix A( fouriertp, 1, ( truncation + 1 ) * 2 );
            eckit::linalg::Matrix B( scl_fourier_tp, ( truncation + 1 ) * 2, nb_fields );
            eckit::linalg::Matrix C( gp_opt, 1, nb_fields );
            linalg_.gemm( A, B, C );
            for ( int j = 0; j < nb_fields; j++ ) {
                gp_fields[ip + j * grid_.size()] = gp_opt[j];
            }
        }
        // Computing u,v from U,V:
        {
            if ( nb_vordiv_fields > 0 ) {
                //ATLAS_TRACE( "u,v from U,V" );
                const double coslat = std::cos( lat );
                for ( int j = 0; j < 2 * nb_vordiv_fields && j < nb_fields; j++ ) {
                    gp_fields[ip + j * grid_.size()] /= coslat;
                }
            }
        }
        ++ip;
    }
    free_aligned( legendre );
    free_aligned( scl_fourier );
    free_aligned( scl_fourier_tp );
    free_aligned( fouriertp );
    free_aligned( gp_opt );
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
void TransLocal::invtrans_uv( const int truncation, const int nb_scalar_fields, const int nb_vordiv_fields,
                              const double scalar_spectra[], double gp_fields[],
                              const eckit::Configuration& config ) const {
    if ( nb_scalar_fields > 0 ) {
        int nb_fields = nb_scalar_fields;

        // Transform
        if ( StructuredGrid( grid_ ) && not grid_.projection() ) {
            auto g = StructuredGrid( grid_ );
            ATLAS_TRACE( "invtrans_uv structured" );
            int nlats            = g.ny();
            int nlons            = g.nxmax();
            int size_fourier_max = nb_fields * 2 * nlats;
            double* scl_fourier;
            alloc_aligned( scl_fourier, size_fourier_max * ( truncation_ + 1 ) );

            // ATLAS-159 workaround begin
            for ( int i = 0; i < size_fourier_max * ( truncation_ + 1 ); ++i ) {
                scl_fourier[i] = 0.;
            }
            // ATLAS-159 workaround end

            // Legendre transformation:
            invtrans_legendre( truncation, nlats, nb_scalar_fields, nb_vordiv_fields, scalar_spectra, scl_fourier,
                               config );

            // Fourier transformation:
            if ( RegularGrid( gridGlobal_ ) ) {
                invtrans_fourier_regular( nlats, nlons, nb_fields, scl_fourier, gp_fields, config );
            }
            else {
                invtrans_fourier_reduced( nlats, g, nb_fields, scl_fourier, gp_fields, config );
            }

            // Computing u,v from U,V:
            {
                if ( nb_vordiv_fields > 0 ) {
                    ATLAS_TRACE( "compute u,v from U,V" );
                    std::vector<double> coslatinvs( nlats );
                    for ( idx_t j = 0; j < nlats; ++j ) {
                        double lat = g.y( j );
                        if ( lat > latPole ) { lat = latPole; }
                        if ( lat < -latPole ) { lat = -latPole; }
                        double coslat = std::cos( lat * util::Constants::degreesToRadians() );
                        coslatinvs[j] = 1. / coslat;
                        //Log::info() << "lat=" << g.y( j ) << " coslat=" << coslat << std::endl;
                    }
                    int idx = 0;
                    for ( idx_t jfld = 0; jfld < 2 * nb_vordiv_fields && jfld < nb_fields; jfld++ ) {
                        for ( idx_t jlat = 0; jlat < g.ny(); jlat++ ) {
                            for ( idx_t jlon = 0; jlon < g.nx( jlat ); jlon++ ) {
                                gp_fields[idx] *= coslatinvs[jlat];
                                idx++;
                            }
                        }
                    }
                }
            }
            free_aligned( scl_fourier );
        }
        else {
            if ( unstruct_precomp_ ) {
                invtrans_unstructured_precomp( truncation, nb_scalar_fields, nb_vordiv_fields, scalar_spectra,
                                               gp_fields, config );
            }
            else {
                invtrans_unstructured( truncation, nb_scalar_fields, nb_vordiv_fields, scalar_spectra, gp_fields,
                                       config );
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans( const int nb_vordiv_fields, const double vorticity_spectra[],
                           const double divergence_spectra[], double gp_fields[],
                           const eckit::Configuration& config ) const {
    invtrans( 0, nullptr, nb_vordiv_fields, vorticity_spectra, divergence_spectra, gp_fields, config );
}

// --------------------------------------------------------------------------------------------------------------------

void extend_truncation( const int old_truncation, const int nb_fields, const double old_spectra[],
                        double new_spectra[] ) {
    int k = 0, k_old = 0;
    for ( int m = 0; m <= old_truncation + 1; m++ ) {             // zonal wavenumber
        for ( int n = m; n <= old_truncation + 1; n++ ) {         // total wavenumber
            for ( int imag = 0; imag < 2; imag++ ) {              // imaginary/real part
                for ( int jfld = 0; jfld < nb_fields; jfld++ ) {  // field
                    if ( m == old_truncation + 1 || n == old_truncation + 1 ) { new_spectra[k++] = 0.; }
                    else {
                        new_spectra[k++] = old_spectra[k_old++];
                    }
                }
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans( const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                           const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                           const eckit::Configuration& config ) const {
    int nb_gp = grid_.size();
    if ( nb_vordiv_fields > 0 ) {
        // collect all spectral data into one array "all_spectra":
        ATLAS_TRACE( "TransLocal::invtrans" );
        int nb_vordiv_spec_ext = 2 * legendre_size( truncation_ + 1 ) * nb_vordiv_fields;
        std::vector<double> U_ext;
        std::vector<double> V_ext;
        std::vector<double> scalar_ext;
        if ( nb_vordiv_fields > 0 ) {
            std::vector<double> vorticity_spectra_extended( nb_vordiv_spec_ext );
            std::vector<double> divergence_spectra_extended( nb_vordiv_spec_ext );
            U_ext.resize( nb_vordiv_spec_ext );
            V_ext.resize( nb_vordiv_spec_ext );

            {
                ATLAS_TRACE( "extend vordiv" );
                // increase truncation in vorticity_spectra and divergence_spectra:
                extend_truncation( truncation_, nb_vordiv_fields, vorticity_spectra,
                                   vorticity_spectra_extended.data() );
                extend_truncation( truncation_, nb_vordiv_fields, divergence_spectra,
                                   divergence_spectra_extended.data() );
            }

            {
                ATLAS_TRACE( "vordiv to UV" );
                // call vd2uv to compute u and v in spectral space
                trans::VorDivToUV vordiv_to_UV_ext( truncation_ + 1, option::type( "local" ) );
                vordiv_to_UV_ext.execute( nb_vordiv_spec_ext, nb_vordiv_fields, vorticity_spectra_extended.data(),
                                          divergence_spectra_extended.data(), U_ext.data(), V_ext.data() );
            }
        }
        if ( nb_scalar_fields > 0 ) {
            int nb_scalar_ext = 2 * legendre_size( truncation_ + 1 ) * nb_scalar_fields;
            scalar_ext.resize( nb_scalar_ext );
            extend_truncation( truncation_, nb_scalar_fields, scalar_spectra, scalar_ext.data() );
        }
        int nb_all_fields = 2 * nb_vordiv_fields + nb_scalar_fields;
        int nb_all_size   = 2 * legendre_size( truncation_ + 1 ) * nb_all_fields;
        std::vector<double> all_spectra( nb_all_size );
        int k = 0, i = 0, j = 0, l = 0;
        {
            ATLAS_TRACE( "merge all spectra" );
            for ( int m = 0; m <= truncation_ + 1; m++ ) {                       // zonal wavenumber
                for ( int n = m; n <= truncation_ + 1; n++ ) {                   // total wavenumber
                    for ( int imag = 0; imag < 2; imag++ ) {                     // imaginary/real part
                        for ( int jfld = 0; jfld < nb_vordiv_fields; jfld++ ) {  // vorticity fields
                            all_spectra[k++] = U_ext[i++];
                        }
                        for ( int jfld = 0; jfld < nb_vordiv_fields; jfld++ ) {  // divergence fields
                            all_spectra[k++] = V_ext[j++];
                        }
                        for ( int jfld = 0; jfld < nb_scalar_fields; jfld++ ) {  // scalar fields
                            all_spectra[k++] = scalar_ext[l++];
                        }
                    }
                }
            }
        }
        int nb_vordiv_size = 2 * legendre_size( truncation_ + 1 ) * nb_vordiv_fields;
        int nb_scalar_size = 2 * legendre_size( truncation_ + 1 ) * nb_scalar_fields;
        ATLAS_ASSERT( k == nb_all_size );
        ATLAS_ASSERT( i == nb_vordiv_size );
        ATLAS_ASSERT( j == nb_vordiv_size );
        ATLAS_ASSERT( l == nb_scalar_size );
        invtrans_uv( truncation_ + 1, nb_all_fields, nb_vordiv_fields, all_spectra.data(), gp_fields, config );
    }
    else {
        if ( nb_scalar_fields > 0 ) {
            invtrans_uv( truncation_, nb_scalar_fields, 0, scalar_spectra, gp_fields + 2 * nb_gp * nb_vordiv_fields,
                         config );
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans( const Field& gpfield, Field& spfield, const eckit::Configuration& config ) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans( const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& config ) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                                       const eckit::Configuration& config ) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                           const eckit::Configuration& ) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                           double divergence_spectra[], const eckit::Configuration& ) const {
    ATLAS_NOTIMPLEMENTED;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
