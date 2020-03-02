/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "eckit/filesystem/PathName.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/DataType.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/detail/PointCloudIO.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace output {
namespace detail {

// ------------------------------------------------------------------

namespace {

std::string sanitize_field_name( const std::string& s ) {
    // replace non-printable characters, then trim right & left
    std::string r( s );
    std::replace_if( r.begin(), r.end(), ::isspace, '_' );
    r.erase( r.find_last_not_of( '_' ) + 1 );
    r.erase( 0, r.find_first_not_of( '_' ) );
    if ( !r.length() ) {
        r = "_";
    }
    return r;
}

}  // end anonymous namespace

// ------------------------------------------------------------------

Mesh PointCloudIO::read( const eckit::PathName& path, std::vector<std::string>& vfnames ) {
    const std::string msg( "PointCloudIO::read: " );

    Log::debug() << "PointCloudIO reading " << path << std::endl;

    Mesh mesh;

    vfnames.clear();

    {
        std::string line;
        std::istringstream iss;
        std::ostringstream oss;
        size_t nb_pts;      // # points/nodes
        size_t nb_columns;  // # data columns (lon,lat,field1,field2,...)
        size_t nb_fld;      // # fields (nb_fld := nb_columns-2)

        // open file and read all of header & data
        std::ifstream f( path.asString().c_str() );
        if ( !f.is_open() ) {
            throw_CantOpenFile( path.asString() );
        }

        // header, part 1:
        // determine number of rows/columns
        // (read all of line, look for "PointCloudIO" signature, nb_pts, nb_columns,
        // ...)
        std::getline( f, line );
        iss.str( line );
        iss >> line >> nb_pts >> nb_columns;
        if ( line != "PointCloudIO" ) {
            std::stringstream errmsg;
            errmsg << msg << "beginning of file `" << path << "` not found (expected: PointCloudIO, got: " << line
                   << ")";
            throw_Exception( errmsg.str(), Here() );
        }
        if ( nb_pts == 0 ) {
            throw_AssertionFailed( msg + " invalid number of points (failed: nb_pts>0)" );
        }
        if ( nb_columns < 2 ) {
            throw_AssertionFailed( msg + " invalid number of columns (failed: nb_columns>=2)" );
        }

        mesh.nodes().resize( static_cast<idx_t>( nb_pts ) );

        mesh::Nodes& nodes                  = mesh.nodes();
        array::ArrayView<double, 2> xy      = array::make_view<double, 2>( nodes.xy() );
        array::ArrayView<double, 2> lonlat  = array::make_view<double, 2>( nodes.lonlat() );
        array::ArrayView<gidx_t, 1> glb_idx = array::make_view<gidx_t, 1>( nodes.global_index() );
        array::ArrayView<int, 1> part       = array::make_view<int, 1>( nodes.partition() );
        part.assign( 0 );
        // header, part 2:
        // determine columns' labels
        // (check end of first line for possible column labels, starting from
        // defaults)

        vfnames.resize( nb_columns );
        for ( size_t j = 0; j < nb_columns; ++j ) {
            std::stringstream name;
            name.str( "column_" );
            name << ( j + 1 );
            vfnames[j] = ( iss && iss >> line ) ? sanitize_field_name( line ) : name.str();
        }

        // (preallocate data, and define fields without the first two columns
        // (lon,lat) because they are treated differently)
        vfnames.erase( vfnames.begin(), vfnames.begin() + 2 );
        nb_fld = nb_columns - 2;  // always >= 0, considering previous check

        std::vector<array::ArrayView<double, 1>> fields;
        for ( size_t j = 0; j < nb_fld; ++j ) {
            fields.emplace_back( array::make_view<double, 1>( nodes.add( Field(
                vfnames[j], array::make_datatype<double>(), array::make_shape( static_cast<idx_t>( nb_pts ) ) ) ) ) );
        }

        size_t i, j;  // (index for node/row and field/column, out of scope to check
                      // at end of loops)
        for ( i = 0; f && i < nb_pts; ++i ) {
            std::getline( f, line );
            iss.clear();
            iss.str( line );

            // NOTE always expects (lon,lat) order, maybe make it configurable?
            PointXY pxy;
            iss >> pxy.x() >> pxy.y();

            xy( i, (size_t)XX )      = pxy.x();
            xy( i, (size_t)YY )      = pxy.y();
            lonlat( i, (size_t)LON ) = pxy.x();
            lonlat( i, (size_t)LAT ) = pxy.y();
            glb_idx( i )             = i + 1;

            for ( j = 0; iss && j < nb_fld; ++j ) {
                iss >> fields[j]( i );
            }
            if ( j < nb_fld ) {
                oss << " Invalid number of fields in data section, on line " << ( i + 1 ) << ", read " << j
                    << " fields, expected " << nb_fld << ".";
                throw_AssertionFailed( msg + oss.str() );
            }
        }
        if ( i < nb_pts ) {
            oss << " Invalid number of lines in data section, read " << ( i ) << " lines, expected " << nb_pts << ".";
            throw_AssertionFailed( msg + oss.str() );
        }

        f.close();
    }

    return mesh;
}

Mesh PointCloudIO::read( const eckit::PathName& path ) {
    std::vector<std::string> vfnames;
    return read( path, vfnames );
}

void PointCloudIO::write( const eckit::PathName& path, const Mesh& mesh ) {
    const std::string msg( "PointCloudIO::write: " );

    Log::debug() << "PointCloudIO writing " << path << std::endl;

    // operate in mesh function space, creating transversing data structures
    // @warning: several copy operations here

    const mesh::Nodes& nodes = mesh.nodes();

    auto lonlat = array::make_view<double, 2>( nodes.lonlat() );
    if ( !lonlat.size() ) {
        throw_Exception( msg + "invalid number of points (failed: nb_pts>0)" );
    }

    // get the fields (sanitized) names and values
    // (bypasses fields ("lonlat"|"lonlat") as shape(1)!=1)
    std::vector<std::string> vfnames;
    std::vector<array::ArrayView<const double, 1>> vfvalues;
    for ( idx_t i = 0; i < nodes.nb_fields(); ++i ) {
        const Field& field = nodes.field( i );
        if ( ( ( field.rank() == 1 && field.shape( 0 ) == lonlat.shape( 0 ) ) ||
               ( field.rank() == 2 && field.shape( 0 ) == lonlat.shape( 0 ) && field.shape( 1 ) == 1 ) ) &&
             field.datatype() == array::DataType::real64() ) {  // FIXME: no support for non-double types!
            vfnames.push_back( sanitize_field_name( field.name() ) );
            vfvalues.push_back( array::make_view<double, 1>( field ) );
        }
    }

    std::ofstream f( path.asString().c_str() );
    if ( !f.is_open() ) {
        throw_CantOpenFile( path.asString() );
    }

    const size_t Npts = lonlat.shape( 0 );
    const size_t Nfld = vfvalues.size();

    // header
    f << "PointCloudIO\t" << Npts << '\t' << ( 2 + Nfld ) << "\tlon\tlat";
    for ( size_t j = 0; j < Nfld; ++j ) {
        f << '\t' << vfnames[j];
    }
    f << '\n';

    // data
    for ( size_t i = 0; i < Npts; ++i ) {
        f << lonlat( i, (size_t)0 ) << '\t' << lonlat( i, (size_t)1 );
        for ( size_t j = 0; j < Nfld; ++j ) {
            f << '\t' << vfvalues[j]( i );
        }
        f << '\n';
    }

    f.close();
}

void PointCloudIO::write( const eckit::PathName& path, const FieldSet& fieldset,
                          const functionspace::NodeColumns& function_space ) {
    const std::string msg( "PointCloudIO::write: " );

    Log::debug() << "PointCloudIO writing " << path << std::endl;

    // operate in field sets with same grid and consistent size(s), creating
    // transversing data structures
    // @warning: several copy operations here

    ATLAS_ASSERT( fieldset.size() );

    array::ArrayView<double, 2> lonlat = array::make_view<double, 2>( function_space.nodes().xy() );
    if ( !lonlat.size() ) {
        throw_Exception( msg + "invalid number of points (failed: nb_pts>0)" );
    }

    // get the fields (sanitized) names and values
    // (bypasses fields ("lonlat"|"lonlat") as shape(1)!=1)
    std::vector<std::string> vfnames;
    std::vector<array::ArrayView<const double, 1>> vfvalues;
    for ( idx_t i = 0; i < fieldset.size(); ++i ) {
        const Field& field = fieldset[i];
        if ( field.shape( 0 ) == lonlat.shape( 0 ) && field.rank() == 1 &&
             field.name() != "glb_idx" )  // FIXME: no support for non-int types!
        {
            vfnames.push_back( sanitize_field_name( field.name() ) );
            vfvalues.push_back( array::make_view<double, 1>( field ) );
        }
    }

    std::ofstream f( path.asString().c_str() );
    if ( !f.is_open() ) {
        throw_CantOpenFile( path.asString() );
    }
    const size_t Npts = lonlat.shape( 0 ), Nfld = vfvalues.size();

    // header
    f << "PointCloudIO\t" << Npts << '\t' << ( 2 + Nfld ) << "\tlon\tlat";
    for ( size_t j = 0; j < Nfld; ++j ) {
        f << '\t' << vfnames[j];
    }
    f << '\n';

    f.precision( std::numeric_limits<double>::digits10 );

    // data
    for ( size_t i = 0; i < Npts; ++i ) {
        f << lonlat( i, (size_t)0 ) << '\t' << lonlat( i, (size_t)1 );
        for ( size_t j = 0; j < Nfld; ++j ) {
            f << '\t' << vfvalues[j]( i );
        }
        f << '\n';
    }

    f.close();
}

void PointCloudIO::write( const eckit::PathName& path, const std::vector<PointLonLat>& pts ) {
    Log::debug() << "PointCloudIO writing " << path << std::endl;

    std::ofstream f( path.asString().c_str() );
    if ( !f.is_open() ) {
        throw_CantOpenFile( path.asString() );
    }

    // header
    f << "PointCloudIO\t" << pts.size() << '\t' << 2 << "\tlon\tlat\n";

    // data
    for ( size_t i = 0; i < pts.size(); ++i ) {
        f << pts[i].lon() << '\t' << pts[i].lat() << '\n';
    }

    f.close();
}

void PointCloudIO::write( const eckit::PathName& path, const std::vector<double>& lon, const std::vector<double>& lat,
                          const std::vector<std::vector<double>*>& vfvalues, const std::vector<std::string>& vfnames ) {
    Log::debug() << "PointCloudIO writing " << path << std::endl;

    const std::string msg( "PointCloudIO::write: " );
    const size_t Npts( lon.size() ), Nfld( vfvalues.size() );
    if ( Npts != lat.size() ) {
        throw_Exception( msg + "number of points inconsistent (failed: #lon == #lat)" );
    }
    if ( Nfld != vfnames.size() ) {
        throw_Exception( msg + "number of fields inconsistent (failed: #vfvalues == #vfnames)" );
    }
    for ( size_t j = 0; j < Nfld; ++j ) {
        if ( Npts != vfvalues[j]->size() ) {
            throw_Exception( msg +
                             "number of points inconsistent (failed: "
                             "#lon == #lat == #*vfvalues[])" );
        }
    }

    std::ofstream f( path.asString().c_str() );
    if ( !f.is_open() ) {
        throw_CantOpenFile( path.asString() );
    }

    // header
    f << "PointCloudIO\t" << Npts << '\t' << ( 2 + Nfld ) << "\tlon\tlat";
    for ( size_t j = 0; j < Nfld; ++j ) {
        f << '\t' << sanitize_field_name( vfnames[j] );
    }
    f << '\n';

    // data
    for ( size_t i = 0; i < Npts; ++i ) {
        f << lon[i] << '\t' << lat[i];
        for ( size_t j = 0; j < Nfld; ++j ) {
            f << '\t' << vfvalues[j]->operator[]( i );
        }
        f << '\n';
    }

    f.close();
}

void PointCloudIO::write( const eckit::PathName& path, const int& nb_pts, const double* lon, const double* lat,
                          const int& nb_fld, const double** afvalues, const char** afnames ) {
    Log::debug() << "PointCloudIO writing " << path << std::endl;


    const std::string msg( "PointCloudIO::write: " );

    const size_t Npts( nb_pts > 0 ? nb_pts : 0 ), Nfld( nb_fld > 0 && afvalues && afnames ? nb_fld : 0 );
    if ( !Npts ) {
        throw_Exception( msg + "invalid number of points (nb_nodes)" );
    }
    if ( !lon ) {
        throw_Exception( msg + "invalid array describing longitude (lon)" );
    }
    if ( !lat ) {
        throw_Exception( msg + "invalid array describing latitude (lat)" );
    }

    std::ofstream f( path.asString().c_str() );
    if ( !f.is_open() ) {
        throw_CantOpenFile( path.asString() );
    }

    // header
    f << "PointCloudIO\t" << Npts << '\t' << ( 2 + Nfld ) << "\tlon\tlat";
    for ( size_t j = 0; j < Nfld; ++j ) {
        f << '\t' << sanitize_field_name( afnames[j] );
    }
    f << '\n';

    // data
    for ( size_t i = 0; i < Npts; ++i ) {
        f << lon[i] << '\t' << lat[i];
        for ( size_t j = 0; j < Nfld; ++j ) {
            f << '\t' << afvalues[j][i];
        }
        f << '\n';
    }

    f.close();
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

// PointCloudIO* atlas__pointcloud__new()
// { return new PointCloudIO(); }
//
//
// void atlas__pointcloud__delete (PointCloudIO* This)
// { delete This; }
//
//
// mesh::Mesh* atlas__pointcloud__read (PointCloudIO* This, char* file_path)
// { return This->read(file_path); }
//
//
// mesh::Mesh* atlas__read_pointcloud (char* file_path)
// { return PointCloudIO::read(file_path); }
//
//
// void atlas__write_pointcloud_fieldset (char* file_path, const
// field::FieldSetImpl* fieldset, const functionspace::detail::NodeColumns*
// functionspace)
// { PointCloudIO::write(file_path, *fieldset, *functionspace); }

// ------------------------------------------------------------------

}  // namespace detail
}  // namespace output
}  // namespace atlas
