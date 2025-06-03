/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#pragma once

#include "StructuredInterpolation2D.h"

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <thread>

#include "eckit/filesystem/LocalPathName.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Point.h"

// With nvidia/22.1, following error is encountered during interpolation.
//     NVC++-S-0000-Internal compiler error. BAD sptr in var_refsym       0
// It is observed to have been fixed with nvidia/22.11
// Disabling OpenMP in this routine for nvidia < 22.11 seems to fix things
#if defined(__NVCOMPILER)
#if (__NVCOMPILER_MAJOR__*100) + __NVCOMPILER_MINOR__ < 2211
#warning "Disabled OpenMP for StructuredInterpolation2D due to internal compiler error"
#undef atlas_omp_parallel
#define atlas_omp_parallel
#undef atlas_omp_for
#define atlas_omp_for for
#endif
#endif



namespace atlas {
namespace interpolation {
namespace method {

namespace structured2d {

namespace {

    inline std::string search_replace(const std::string& in, const std::string& search, const std::string& replace) {
        std::string out = in;
        int pos = out.find(search);
        while (pos != std::string::npos) {
            out.erase(pos, search.length());
            out.insert(pos, replace);
            pos = out.find(search, pos + replace.length());
        }
        return out;
    }

    inline std::string left_padded_str(int x, int max = 0) {
        auto digits = [](long x) -> long { return std::floor(std::log10(std::max(1l, x))) + 1l; };
        if (max) {
            std::ostringstream ss;
            ss << std::setw(digits(max)) << std::setfill('0') << x;
            return ss.str();
        }
        return std::to_string(x);
    }

    inline std::string filename( const std::string& path ){
        return search_replace(path, "%p", left_padded_str(mpi::rank(),mpi::size()));
    }

    template <typename LonLat>
    inline std::string to_str (const std::vector<idx_t>& points, LonLat lonlat) {
        std::ostringstream out;
        out << "[\n";
        for( size_t j=0; j<points.size(); ++j ) {
            PointLonLat p = lonlat(points[j]);
            out << "  [" << p.lon() << "," << p.lat() << "]";
            if( j < points.size() - 1 ) {
                out << ",\n";
            }
        }
        out << "\n]";
        return out.str();
    };

    inline void output_source_partition_polygon(const std::string& path, functionspace::StructuredColumns src_fs, idx_t halo) {
        auto& polygon = src_fs.polygon(halo);
        std::ofstream f(filename(path));
        f << "[\n  " << polygon.json() << "\n]";
    }

    inline void output_all_source_partition_polygons(const std::string& path, functionspace::StructuredColumns src_fs, idx_t halo) {
        auto& polygon = src_fs.polygon(halo);
        util::PartitionPolygons polygons;
        polygon.allGather(polygons);
        if( mpi::rank() == 0 ) {
            std::ofstream f(path);
            f << polygons.json();
        }
    }

    inline void output_source_grid_points(const std::string& path, functionspace::StructuredColumns src_fs, idx_t halo) {
        std::ofstream f(filename(path));
        f << "[\n";
        if( halo == 0 ) {
            idx_t size = src_fs.sizeOwned();
            idx_t n = 0;
            for( idx_t j=src_fs.j_begin(); j<src_fs.j_end(); ++j ) {
                for( idx_t i=src_fs.i_begin(j); i<src_fs.i_end(j); ++i ) {
                    auto p = src_fs.compute_xy(i,j);
                    f << "  [" << p.x() << "," << p.y() << "]";
                    if( (++n) != size ) {
                        f << ",\n";
                    }
                }
            }
        }
        else {
            idx_t size = src_fs.sizeHalo();
            idx_t n = 0;
            for( idx_t j=src_fs.j_begin_halo(); j<src_fs.j_end_halo(); ++j ) {
                for( idx_t i=src_fs.i_begin_halo(j); i<src_fs.i_end_halo(j); ++i ) {
                    auto p = src_fs.compute_xy(i,j);
                    f << "  [" << p.x() << "," << p.y() << "]";
                    if( (++n) != size ) {
                        f << ",\n";
                    }
                }
            }
        }
        f << "\n]";
    }


    template <typename LonLat>
    inline void output_target_points(const std::string& path, const std::vector<idx_t>& points, LonLat lonlat) {
        std::ofstream f(filename(path));
        f << to_str(points,lonlat);
    }

    template <typename LonLat, typename Interpolation>
    inline void handle_failed_points(const Interpolation& interpolation, const std::vector<idx_t>& failed_points, LonLat lonlat ) {
        const auto src_fs = functionspace::StructuredColumns( interpolation.source() );
        size_t num_failed_points{0};
        mpi::comm().allReduce(failed_points.size(), num_failed_points, eckit::mpi::sum());
        if( num_failed_points > 0 ) {
            std::string halo_str = std::to_string(src_fs.halo());
            output_source_partition_polygon("atlas_source_partition_polygons_halo_0_p%p.json",src_fs,0);
            output_source_partition_polygon("atlas_source_partition_polygons_halo_" + halo_str + "_p%p.json",src_fs,src_fs.halo());
            output_target_points("atlas_target_failed_points_p%p.json", failed_points, lonlat);
            idx_t my_rank = mpi::rank();
            for( idx_t p=0; p<mpi::size(); ++p ) {
                if( p == my_rank && failed_points.size() ) {
                    Log::error() << "Failed to interpolate " << failed_points.size() << " points on rank " << p
                                 << ".\nSee " << filename("atlas_target_failed_points_p%p.json") ;
                    if( failed_points.size() <= 20 ) {
                        Log::error() << " :\n" << to_str(failed_points, lonlat);
                    }
                    Log::error() << std::endl;
                }
                mpi::comm().barrier();
            }
            output_all_source_partition_polygons("atlas_source_partition_polygons_halo_0.json",src_fs,0);
            output_all_source_partition_polygons("atlas_source_partition_polygons_halo_"+halo_str+".json",src_fs,src_fs.halo());
            // output_source_grid_points("atlas_source_grid_points_halo_0_p%p.json",src_fs,0);
            // output_source_grid_points("atlas_source_grid_points_halo_"+halo_str+"_p%p.json",src_fs,src_fs.halo());
            {
                std::ofstream f("atlas_interpolation_info.json");
                util::Config src_fs_config;
                src_fs_config.set("source.functionspace","StructuredColumns");
                src_fs_config.set("source.grid",src_fs.grid().spec());
                src_fs_config.set("source.halo",src_fs.halo());
                src_fs_config.set("source.distribution",src_fs.distribution());
                src_fs_config.set("source.partitions",src_fs.nb_parts());
                src_fs_config.set("interpolation","StructuredInterpolation2D<"+interpolation.kernel().className()+">");
                f << src_fs_config.json();
            }
            std::string p_range = "{"+left_padded_str(0,mpi::size()-1)+".."+left_padded_str(mpi::size()-1,mpi::size()-1)+"}";
            Log::info() << "Dumped files to " << eckit::LocalPathName::cwd() << ": \n"
                        << "    atlas_target_failed_points_p" << p_range << ".json\n"
                        << "    atlas_source_partition_polygons_halo_0_p" << p_range << ".json\n"
                        << "    atlas_source_partition_polygons_halo_0.json\n"
                        << "    atlas_source_partition_polygons_halo_" << src_fs.halo() << "_p" << p_range << ".json\n"
                        << "    atlas_source_partition_polygons_halo_" << src_fs.halo() << ".json\n"
                        << "    atlas_interpolation_info.json\n"
                        << std::endl;

            std::ostringstream err;
            err << "StructuredInterpolation2D<" << interpolation.kernel().className() << "> failed for "
                << num_failed_points << " points with source halo="<<src_fs.halo()
                << ". Try increasing the source halo. Files have been written for debugging purpose to ["
                << eckit::LocalPathName::cwd() << "].";
            ATLAS_THROW_EXCEPTION(err.str());
        }
    }
}
}


template <typename Kernel>
double StructuredInterpolation2D<Kernel>::convert_units_multiplier( const Field& field ) {
    std::string units = field.metadata().getString( "units", "degrees" );
    if ( units == "degrees" ) {
        return 1.;
    }
    if ( units == "radians" ) {
        return 180. / M_PI;
    }
    ATLAS_NOTIMPLEMENTED;
}

template <typename Kernel>
StructuredInterpolation2D<Kernel>::StructuredInterpolation2D( const Method::Config& config ) :
    Method( config ),
    verbose_{false},
    limiter_{false},
    matrix_free_{false} {
    config.get( "verbose", verbose_ );
    config.get( "limiter", limiter_ );
    config.get( "matrix_free", matrix_free_ );
    if (limiter_ && not matrix_free_) {
        ATLAS_THROW_EXCEPTION("Cannot apply configuration 'limiter=true' and 'matrix_free=false' together");
    }
}


template <typename Kernel>
void StructuredInterpolation2D<Kernel>::do_setup( const Grid& source, const Grid& target, const Cache& ) {
    ATLAS_TRACE( "StructuredInterpolation2D<" + Kernel::className() + ">::do_setup(Grid source, Grid target)" );
    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }


    ATLAS_ASSERT( StructuredGrid( source ) );
    FunctionSpace source_fs =
        functionspace::StructuredColumns( source, option::halo( std::max<idx_t>( kernel_->stencil_halo(), 1 ) ) );
    // guarantee "1" halo for pole treatment!
    FunctionSpace target_fs = functionspace::PointCloud( target );

    do_setup( source_fs, target_fs );
}


template <typename Kernel>
void StructuredInterpolation2D<Kernel>::do_setup( const FunctionSpace& source, const FunctionSpace& target, const Cache& cache) {
    ATLAS_TRACE( "StructuredInterpolation2D<" + Kernel::className() + ">::do_setup(FunctionSpace source, FunctionSpace target)" );
    if (! matrix_free_ && interpolation::MatrixCache(cache)) {
        setMatrix(cache);
        source_ = source;
        target_ = target;
        ATLAS_ASSERT(matrix().rows() == target.size());
        ATLAS_ASSERT(matrix().cols() == source.size());
        return;
    }
    else {
        do_setup(source, target);
    }
}

template <typename Kernel>
void StructuredInterpolation2D<Kernel>::do_setup( const FunctionSpace& source, const FunctionSpace& target ) {
    ATLAS_TRACE( "StructuredInterpolation2D<" + Kernel::className() + ">::do_setup(FS source, FS target)" );

    source_ = source;
    target_ = target;

    if ( functionspace::NodeColumns tgt = target ) {
        target_lonlat_ = tgt.mesh().nodes().lonlat();
        target_ghost_  = tgt.mesh().nodes().ghost();
    }
    else if ( functionspace::PointCloud tgt = target ) {
        target_lonlat_ = tgt.lonlat();
        target_ghost_  = tgt.ghost();
    }
    else if ( functionspace::StructuredColumns tgt = target ) {
        target_lonlat_ = tgt.xy();
        target_ghost_  = tgt.ghost();
    }
    else {
        throw_NotImplemented(
            "Only interpolation to functionspaces NodeColumns, PointCloud or StructuredColumns are implemented",
            Here() );
    }

    setup( source );
}

template <typename Kernel>
void StructuredInterpolation2D<Kernel>::do_setup( const FunctionSpace& source, const Field& target ) {
    ATLAS_TRACE( "StructuredInterpolation<" + Kernel::className() + ">::do_setup(FunctionSpace source, Field target)" );

    source_ = source;

    if ( target.functionspace() ) {
        target_ = target.functionspace();
    }

    target_lonlat_ = target;

    setup( source );
}

template <typename Kernel>
void StructuredInterpolation2D<Kernel>::do_setup( const FunctionSpace& source, const FieldSet& target ) {
    ATLAS_TRACE( "StructuredInterpolation<" + Kernel::className() + ">::do_setup(FunctionSpace source,FieldSet target)" );

    source_ = source;

    ATLAS_ASSERT( target.size() >= 2 );
    if ( target[0].functionspace() ) {
        target_ = target[0].functionspace();
    }

    target_lonlat_fields_ = target;

    setup( source );
}

template <typename Kernel>
void StructuredInterpolation2D<Kernel>::print( std::ostream& ) const {
    ATLAS_NOTIMPLEMENTED;
}


template <typename Kernel>
void StructuredInterpolation2D<Kernel>::setup( const FunctionSpace& source ) {
    using namespace structured2d;

    kernel_.reset( new Kernel( source, util::Config( "limiter", limiter_ ) ) );

    if ( functionspace::StructuredColumns( source ).halo() < 1 ) {
        throw_Exception( "The source functionspace must have (halo >= 1) for pole treatment" );
    }

    out_npts_ = 0;
    if( target_lonlat_ ) {
        convert_units_ = convert_units_multiplier(target_lonlat_);
        out_npts_ = target_lonlat_.shape( 0 );
    }
    else if( target_lonlat_fields_.empty() ) {
        convert_units_ = convert_units_multiplier(target_lonlat_fields_[LON]);
        out_npts_ = target_lonlat_fields_[0].shape(0);
    }

    if ( not matrix_free_ ) {
        ATLAS_TRACE( "Precomputing interpolation matrix" );

        std::vector<idx_t> failed_points;

        auto triplets = kernel_->allocate_triplets( out_npts_ );

        using WorkSpace = typename Kernel::WorkSpace;
        auto interpolate_point = [&]( idx_t n, PointLonLat&& p, WorkSpace& workspace ) -> int {
            try {
                kernel_->insert_triplets( n, p, triplets, workspace );
                return 0;
            }
            catch(const eckit::Exception& e) {}
            if (verbose_) {
                Log::error() << "Could not interpolate point " << n << " :\t" << p << std::endl;
            }
            return 1;
        };

        auto interpolate_omp = [&failed_points,interpolate_point]( idx_t out_npts, auto lonlat, auto ghost) {
            atlas_omp_parallel {
                WorkSpace workspace;
                atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                    if( not ghost(n) ) {
                        if (interpolate_point(n, lonlat(n), workspace) != 0) {
                            atlas_omp_critical {
                                failed_points.emplace_back(n);
                            }
                        }
                    }
                }
            }
        };

        if ( target_lonlat_ ) {
            auto lonlat_view    = array::make_view<double, 2>( target_lonlat_ );
            auto lonlat = [lonlat_view, convert_units = convert_units_] (idx_t n) {
                return PointLonLat{lonlat_view(n,LON) * convert_units, lonlat_view(n,LAT) * convert_units};
            };

            if( out_npts_ != 0 ) {
                if ( target_ghost_ ) {
                    auto ghost     = array::make_view<int, 1>( target_ghost_ );
                    interpolate_omp(out_npts_, lonlat, ghost);
                }
                else {
                    auto no_ghost = [](idx_t n) { return false; };
                    interpolate_omp(out_npts_, lonlat, no_ghost);
                }
            }
            handle_failed_points(*this, failed_points, lonlat);
        }
        else if ( not target_lonlat_fields_.empty() ) {
            const auto lon = array::make_view<double, 1>( target_lonlat_fields_[LON] );
            const auto lat = array::make_view<double, 1>( target_lonlat_fields_[LAT] );
            auto lonlat = [lon, lat, convert_units = convert_units_] (idx_t n) {
                return PointLonLat{lon(n) * convert_units, lat(n) * convert_units};
            };

            if( out_npts_ != 0) {
                if ( target_ghost_ ) {
                    auto ghost     = array::make_view<int, 1>( target_ghost_ );
                    interpolate_omp(out_npts_, lonlat, ghost);
                }
                else {
                    auto no_ghost = [](idx_t n) { return false; };
                    interpolate_omp(out_npts_, lonlat, no_ghost);
                }
            }
            handle_failed_points(*this, failed_points, lonlat);
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }

        // fill sparse matrix
        if( failed_points.empty() && out_npts_) {
            idx_t inp_npts = source.size();
            setMatrix(out_npts_, inp_npts, triplets);
        }
    }
}


template <typename Kernel>
void StructuredInterpolation2D<Kernel>::do_execute( const Field& src_field, Field& tgt_field, Metadata& metadata ) const {
    FieldSet tgt( tgt_field );
    do_execute( FieldSet( src_field ), tgt, metadata );
}


template <typename Kernel>
void StructuredInterpolation2D<Kernel>::do_execute( const FieldSet& src_fields, FieldSet& tgt_fields, Metadata& metadata ) const {
    if ( not matrix_free_ ) {
        Method::do_execute( src_fields, tgt_fields, metadata );
        return;
    }

    ATLAS_TRACE( "StructuredInterpolation<" + Kernel::className() + ">::do_execute()" );

    const idx_t N = src_fields.size();
    ATLAS_ASSERT( N == tgt_fields.size() );

    if ( N == 0 )
        return;

    haloExchange( src_fields );

    array::DataType datatype = src_fields[0].datatype();
    int rank                 = src_fields[0].rank();

    for ( idx_t i = 0; i < N; ++i ) {
        ATLAS_ASSERT( src_fields[i].datatype() == datatype );
        ATLAS_ASSERT( src_fields[i].rank() == rank );
        ATLAS_ASSERT( tgt_fields[i].datatype() == datatype );
        ATLAS_ASSERT( tgt_fields[i].rank() == rank );
    }

    if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 1 ) {
        execute_impl<double, 1>( *kernel_, src_fields, tgt_fields );
    }
    else if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 1 ) {
        execute_impl<float, 1>( *kernel_, src_fields, tgt_fields );
    }
    else if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 2 ) {
        execute_impl<double, 2>( *kernel_, src_fields, tgt_fields );
    }
    else if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 2 ) {
        execute_impl<float, 2>( *kernel_, src_fields, tgt_fields );
    }
    else if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 3 ) {
        execute_impl<double, 3>( *kernel_, src_fields, tgt_fields );
    }
    else if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 3 ) {
        execute_impl<float, 3>( *kernel_, src_fields, tgt_fields );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }

    tgt_fields.set_dirty();
}


template <typename Kernel>
template <typename Value, int Rank>
void StructuredInterpolation2D<Kernel>::execute_impl( const Kernel& kernel, const FieldSet& src_fields,
                                                      FieldSet& tgt_fields ) const {
    using namespace structured2d;

    const idx_t N = src_fields.size();

    std::vector<array::ArrayView<const Value, Rank> > src_view;
    std::vector<array::ArrayView<Value, Rank> > tgt_view;
    src_view.reserve( N );
    tgt_view.reserve( N );

    for ( idx_t i = 0; i < N; ++i ) {
        src_view.emplace_back( array::make_view<Value, Rank>( src_fields[i] ) );
        tgt_view.emplace_back( array::make_view<Value, Rank>( tgt_fields[i] ) );
    }

    using WorkSpace = typename Kernel::WorkSpace;

    auto interpolate_point = [&]( idx_t n, PointLonLat&& p, WorkSpace& workspace ) -> int {
        try {
            kernel.compute_stencil( p.lon(), p.lat(), workspace.stencil );
            kernel.compute_weights( p.lon(), p.lat(), workspace.stencil, workspace.weights );
            kernel.make_valid_stencil( p.lon(), p.lat(), workspace.stencil );
            for ( idx_t i = 0; i < N; ++i ) {
                kernel.interpolate( workspace.stencil, workspace.weights, src_view[i], tgt_view[i], n );
            }
            return 0;
        }
        catch(const eckit::Exception& e) {}
        if (verbose_) {
            Log::error() << "Could not interpolate point " << n << " :\t" << p << std::endl;
        }
        return 1;
    };

    std::vector<idx_t> failed_points;

    auto interpolate_omp = [&failed_points,interpolate_point]( idx_t out_npts, auto lonlat, auto ghost) {
        atlas_omp_parallel {
            WorkSpace workspace;
            atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                if( not ghost(n) ) {
                    if (interpolate_point(n, lonlat(n), workspace) != 0) {
                        atlas_omp_critical {
                            failed_points.emplace_back(n);
                        }
                    }
                }
            }
        }
    };

    if ( target_lonlat_ ) {
        auto lonlat_view    = array::make_view<double, 2>( target_lonlat_ );
        auto lonlat = [lonlat_view, convert_units = convert_units_] (idx_t n) {
            return PointLonLat{lonlat_view(n,LON) * convert_units, lonlat_view(n,LAT) * convert_units};
        };

        if( out_npts_ != 0 ) {
            if ( target_ghost_ ) {
                auto ghost     = array::make_view<int, 1>( target_ghost_ );
                interpolate_omp(out_npts_, lonlat, ghost);
            }
            else {
                auto no_ghost = [](idx_t n) { return false; };
                interpolate_omp(out_npts_, lonlat, no_ghost);
            }
        }
        handle_failed_points(*this, failed_points, lonlat);
    }
    else if ( not target_lonlat_fields_.empty() ) {
        const auto lon = array::make_view<double, 1>( target_lonlat_fields_[LON] );
        const auto lat = array::make_view<double, 1>( target_lonlat_fields_[LAT] );
        auto lonlat = [lon, lat, convert_units = convert_units_] (idx_t n) {
            return PointLonLat{lon(n) * convert_units, lat(n) * convert_units};
        };

        if( out_npts_ != 0) {
            if ( target_ghost_ ) {
                auto ghost     = array::make_view<int, 1>( target_ghost_ );
                interpolate_omp(out_npts_, lonlat, ghost);
            }
            else {
                auto no_ghost = [](idx_t n) { return false; };
                interpolate_omp(out_npts_, lonlat, no_ghost);
            }
        }
        handle_failed_points(*this, failed_points, lonlat);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
