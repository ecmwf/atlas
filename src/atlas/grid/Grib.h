/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_GribWrite_h
#define atlas_GribWrite_h

#include "eckit/grib/GribHandle.h"

#include "atlas/grid/FieldSet.h"
#include "atlas/grid/GridSpec.h"

//------------------------------------------------------------------------------------------------------

namespace eckit { class PathName; }
namespace eckit { class DataHandle; }

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

class Grib {

public: // methods

	static atlas::grid::Grid::Ptr create_grid( eckit::grib::GribHandle& );

	/// Given a Grid, this function will find the closest matching GRIB samples file.
	/// The cloned/new handle of the GRIB sample file is returned.
	/// If no match found a NULL handle is returned.
	static eckit::grib::GribHandle::Ptr create_handle( const grid::Grid&, long edition = 0 );

	/// Given a GridSpec return closest grib samples file.
	/// If no match found returns an empty string
	static std::string grib_sample_file( const grid::GridSpec&, long editionNumber );

	/// Helper function, used locally and in testing
	static void determine_grib_samples_dir(std::vector<std::string>& sample_paths);

	static void write( const atlas::grid::FieldSet& fset, const eckit::PathName& opath  );

	static void write( const atlas::grid::FieldHandle& fh, eckit::DataHandle& dh );

	static eckit::grib::GribHandle::Ptr write( const atlas::grid::FieldHandle& fh );

	static void clone( const atlas::grid::FieldSet& field, const eckit::PathName& src, const eckit::PathName& opath  );

private: // methods

	static void write( const atlas::grid::FieldHandle& field, const eckit::PathName& opath  );

	static void clone( const atlas::grid::FieldHandle& field, const eckit::PathName& gridsec, eckit::DataHandle& );

	static eckit::grib::GribHandle::Ptr clone(const grid::FieldHandle &field, eckit::grib::GribHandle& gridsec );

	/// @todo this function is temporary, until we make an abstract interface to output to different formats
	///       we must learn more about NetCDF, etc...
	///
	static void write_gridspec_to_grib( const GridSpec&, eckit::grib::GribHandle& );

};

//---------------------------------------------------------------------------------------------------------


} // namespace grid
} // namespace atlas

#endif

