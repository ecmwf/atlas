/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// file deepcode ignore MissingOpenCheckOnFile: False positive

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "eckit/filesystem/PathName.h"

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

using atlas::functionspace::NodeColumns;
using atlas::util::Metadata;
using atlas::util::Topology;
using eckit::PathName;

namespace atlas {
namespace output {
namespace detail {

namespace {

class GmshFile : public std::ofstream {
public:
    GmshFile(const PathName& file_path, std::ios_base::openmode mode, int part = static_cast<int>(mpi::rank())) {
        PathName par_path(file_path);
        int mpi_size = static_cast<int>(mpi::size());
        if (mpi::size() == 1 || part == -1) {
            std::ofstream::open(par_path.localPath(), mode);
        }
        else {
            if (mpi::rank() == 0) {
                PathName par_path(file_path);
                std::ofstream par_file(par_path.localPath(), std::ios_base::out);
                for (int p = 0; p < mpi_size; ++p) {
                    PathName loc_path(file_path);
                    // loc_path = loc_path.baseName(false) + "_p" + to_str(p) + ".msh";
                    loc_path = loc_path.baseName(false) + ".msh.p" + std::to_string(p);
                    par_file << "Merge \"" << loc_path << "\";" << std::endl;
                }
                par_file.close();
            }
            PathName path(file_path);
            // path = path.dirName() + "/" + path.baseName(false) + "_p" +
            // to_str(part) + ".msh";
            path = path.dirName() + "/" + path.baseName(false) + ".msh.p" + std::to_string(part);
            std::ofstream::open(path.localPath(), mode);
        }
    }
};

enum GmshElementTypes
{
    LINE  = 1,
    TRIAG = 2,
    QUAD  = 3,
    POINT = 15
};

// ----------------------------------------------------------------------------
void write_header_ascii(std::ostream& out) {
    out << "$MeshFormat\n";
    out << "2.2 0 " << sizeof(double) << "\n";
    out << "$EndMeshFormat\n";
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void write_header_binary(std::ostream& out) {
    out << "$MeshFormat\n";
    out << "2.2 1 " << sizeof(double) << "\n";
    int one = 1;
    out.write(reinterpret_cast<const char*>(&one), sizeof(int));
    out << "\n$EndMeshFormat\n";
}
// ----------------------------------------------------------------------------

namespace {  // anonymous

template <typename T>
array::LocalView<const T, 2> make_level_view(const Field& field, int ndata, int jlev) {
    using namespace array;
    if (field.levels()) {
        if (field.variables()) {
            return make_view<const T, 3>(field).slice(Range::to(ndata), jlev, Range::all());
        }
        else {
            return make_view<const T, 2>(field).slice(Range::to(ndata), jlev, Range::dummy());
        }
    }
    else {
        if (field.variables()) {
            return make_view<const T, 2>(field).slice(Range::to(ndata), Range::all());
        }
        else {
            return make_view<const T, 1>(field).slice(Range::to(ndata), Range::dummy());
        }
    }
}

template <typename Value, typename GlobalIndex, typename IncludeIndex>
void write_level(std::ostream& out, GlobalIndex gidx, const array::LocalView<Value, 2>& data, IncludeIndex include) {
    using value_type = typename std::remove_const<Value>::type;
    int ndata        = data.shape(0);
    int nvars        = data.shape(1);
    if (nvars == 1) {
        for (idx_t n = 0; n < ndata; ++n) {
            if (include(n)) {
                out << gidx(n) << " " << data(n, 0) << "\n";
            }
        }
    }
    else if (nvars <= 3) {
        std::array<value_type, 3> data_vec;
        for (idx_t n = 0; n < ndata; ++n) {
            if (include(n)) {
                for (idx_t v = 0; v < nvars; ++v) {
                    data_vec[v] = data(n, v);
                }
                out << gidx(n);
                for (int v = 0; v < 3; ++v) {
                    out << " " << data_vec[v];
                }
                out << "\n";
            }
        }
    }
    else if (nvars <= 9) {
        std::array<value_type, 9> data_vec;
        if (nvars == 4) {
            for (int n = 0; n < ndata; ++n) {
                if (include(n)) {
                    for (int i = 0; i < 2; ++i) {
                        for (int j = 0; j < 2; ++j) {
                            data_vec[i * 3 + j] = data(n, i * 2 + j);
                        }
                    }
                    out << gidx(n);
                    for (int v = 0; v < 9; ++v) {
                        out << " " << data_vec[v];
                    }
                    out << "\n";
                }
            }
        }
        else if (nvars == 9) {
            for (int n = 0; n < ndata; ++n) {
                if (include(n)) {
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            data_vec[i * 3 + j] = data(n, i * 2 + j);
                        }
                    }
                    out << gidx(n);
                    for (int v = 0; v < 9; ++v) {
                        out << " " << data_vec[v];
                    }
                    out << "\n";
                }
            }
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}
template <typename Value, typename GlobalIndex>
void write_level(std::ostream& out, GlobalIndex gidx, const array::LocalView<Value, 2>& data) {
    write_level(out, gidx, data, [](idx_t) { return true; });
}

std::vector<int> get_levels(int nlev, const Metadata& gmsh_options) {
    std::vector<int> lev;
    std::vector<int> gmsh_levels;
    gmsh_options.get("levels", gmsh_levels);
    if (gmsh_levels.empty() || nlev == 1) {
        lev.resize(nlev);
        for (int ilev = 0; ilev < nlev; ++ilev) {
            lev[ilev] = ilev;
        }
    }
    else {
        lev = gmsh_levels;
    }
    return lev;
}

std::string field_lev(const Field& field, int jlev) {
    if (field.levels()) {
        char str[6] = {0, 0, 0, 0, 0, 0};
        std::sprintf(str, "[%03d]", jlev);
        return std::string(str);
    }
    else {
        return std::string();
    }
}

double field_time(const Field& field) {
    return field.metadata().has("time") ? field.metadata().get<double>("time") : 0.;
}

int field_step(const Field& field) {
    return field.metadata().has("step") ? field.metadata().get<size_t>("step") : 0;
}

int field_vars(int nvars) {
    if (nvars == 1) {
        return nvars;
    }
    else if (nvars <= 3) {
        return 3;
    }
    else if (nvars <= 9) {
        return 9;
    }
    else {
        return nvars;
    }
}

}  // namespace

// ----------------------------------------------------------------------------
template <typename Value>
void write_field_nodes(const Metadata& gmsh_options, const functionspace::NodeColumns& function_space,
                       const Field& field, std::ostream& out) {
    Log::debug() << "writing NodeColumns field " << field.name() << " defined in NodeColumns..." << std::endl;

    bool gather(gmsh_options.get<bool>("gather") && mpi::size() > 1);
    // unused: bool binary( !gmsh_options.get<bool>( "ascii" ) );
    idx_t nlev  = std::max<idx_t>(1, field.levels());
    idx_t ndata = std::min<idx_t>(function_space.nb_nodes(), field.shape(0));
    idx_t nvars = std::max<idx_t>(1, field.variables());
    auto gidx   = array::make_view<gidx_t, 1>(function_space.nodes().global_index());
    Field gidx_glb;
    Field field_glb;
    if (gather) {
        gidx_glb =
            function_space.createField<gidx_t>(option::name("gidx_glb") | option::levels(false) | option::global());
        function_space.gather(function_space.nodes().global_index(), gidx_glb);
        gidx = array::make_view<gidx_t, 1>(gidx_glb);

        field_glb = function_space.createField(field, option::global());
        function_space.gather(field, field_glb);
        ndata = std::min<idx_t>(function_space.nb_nodes_global(), field_glb.shape(0));
    }
    auto missing = field::MissingValue(field);

    std::vector<int> lev = get_levels(nlev, gmsh_options);
    for (size_t ilev = 0; ilev < lev.size(); ++ilev) {
        int jlev = lev[ilev];
        if ((gather && mpi::rank() == 0) || !gather) {
            auto data =
                gather ? make_level_view<Value>(field_glb, ndata, jlev) : make_level_view<Value>(field, ndata, jlev);
            auto include_idx = [&](idx_t n) {
                for (idx_t v = 0; v < nvars; ++v) {
                    if (missing(data(n, v))) {
                        return false;
                    }
                }
                return true;
            };
            idx_t ndata_nonmissing = [&] {
                if (missing) {
                    idx_t c = 0;
                    for (idx_t n = 0; n < ndata; ++n) {
                        c += include_idx(n);
                    }
                    return c;
                }
                return ndata;
            }();

            out << "$NodeData\n";
            out << "1\n";
            out << "\"" << field.name() << field_lev(field, jlev) << "\"\n";
            out << "1\n";
            out << field_time(field) << "\n";
            out << "4\n";
            out << field_step(field) << "\n";
            out << field_vars(nvars) << "\n";
            out << ndata_nonmissing << "\n";
            out << mpi::rank() << "\n";
            if (missing) {
                write_level(out, gidx, data, include_idx);
            }
            else {
                write_level(out, gidx, data);
            }
            out << "$EndNodeData\n";
        }
    }
}

// ----------------------------------------------------------------------------
template <typename Value>
void write_field_nodes(const Metadata& gmsh_options, const functionspace::NoFunctionSpace& function_space,
                       const Field& field, std::ostream& out) {
    Log::debug() << "writing field " << field.name() << " defined without functionspace..." << std::endl;

    // unused: bool binary( !gmsh_options.get<bool>( "ascii" ) );
    idx_t nlev  = std::max<idx_t>(1, field.levels());
    idx_t ndata = field.shape(0);
    idx_t nvars = std::max<idx_t>(1, field.variables());
    auto gidx   = [](idx_t inode) { return inode + 1; };

    auto missing = field::MissingValue(field);

    std::vector<int> lev = get_levels(nlev, gmsh_options);
    for (size_t ilev = 0; ilev < lev.size(); ++ilev) {
        int jlev = lev[ilev];

        auto data        = make_level_view<Value>(field, ndata, jlev);
        auto include_idx = [&](idx_t n) {
            for (idx_t v = 0; v < nvars; ++v) {
                if (missing(data(n, v))) {
                    return false;
                }
            }
            return true;
        };
        idx_t ndata_nonmissing = [&] {
            if (missing) {
                idx_t c = 0;
                for (idx_t n = 0; n < ndata; ++n) {
                    c += include_idx(n);
                }
                return c;
            }
            return ndata;
        }();

        out << "$NodeData\n";
        out << "1\n";
        out << "\"" << field.name() << field_lev(field, jlev) << "\"\n";
        out << "1\n";
        out << field_time(field) << "\n";
        out << "4\n";
        out << field_step(field) << "\n";
        out << field_vars(nvars) << "\n";
        out << ndata_nonmissing << "\n";
        out << mpi::rank() << "\n";
        if (missing) {
            write_level(out, gidx, data, include_idx);
        }
        else {
            write_level(out, gidx, data);
        }
        out << "$EndNodeData\n";
    }
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

void print_field_lev(char field_lev[], int jlev) {
    std::sprintf(field_lev, "[%03d]", jlev);
}

/* unused
void print_field_lev( char field_lev[], long jlev ) {
    std::sprintf( field_lev, "[%03ld]", jlev );
}

void print_field_lev( char field_lev[], unsigned long jlev ) {
    std::sprintf( field_lev, "[%03lu]", jlev );
}
**/

// ----------------------------------------------------------------------------
template <typename DATATYPE>
void write_field_nodes(const Metadata& gmsh_options, const functionspace::StructuredColumns& function_space,
                       const Field& field, std::ostream& out) {
    Log::debug() << "writing StructuredColumns field " << field.name() << "..." << std::endl;

    bool gather(gmsh_options.get<bool>("gather") && mpi::size() > 1);
    // unused: bool binary( !gmsh_options.get<bool>( "ascii" ) );
    idx_t nlev  = std::max<idx_t>(1, field.levels());
    idx_t ndata = std::min<idx_t>(function_space.sizeOwned(), field.shape(0));
    idx_t nvars = std::max<idx_t>(1, field.variables());
    auto gidx   = array::make_view<gidx_t, 1>(function_space.global_index());
    Field gidx_glb;
    Field field_glb;
    if (gather) {
        gidx_glb =
            function_space.createField(function_space.global_index(), option::name("gidx_glb") | option::global());
        function_space.gather(function_space.global_index(), gidx_glb);
        gidx = array::make_view<gidx_t, 1>(gidx_glb);

        field_glb = function_space.createField(field, option::global());
        function_space.gather(field, field_glb);
        ndata = field_glb.shape(0);
    }

    std::vector<int> lev = get_levels(nlev, gmsh_options);
    for (size_t ilev = 0; ilev < lev.size(); ++ilev) {
        int jlev          = lev[ilev];
        char field_lev[6] = {0, 0, 0, 0, 0, 0};

        if (field.levels()) {
            print_field_lev(field_lev, jlev);
        }

        out << "$NodeData\n";
        out << "1\n";
        out << "\"" << field.name() << field_lev << "\"\n";
        out << "1\n";
        out << field_time(field) << "\n";
        out << "4\n";
        out << field_step(field) << "\n";
        out << field_vars(nvars) << "\n";
        out << ndata << "\n";
        out << mpi::rank() << "\n";
        auto data =
            gather ? make_level_view<DATATYPE>(field_glb, ndata, jlev) : make_level_view<DATATYPE>(field, ndata, jlev);
        write_level(out, gidx, data);
        out << "$EndNodeData\n";
    }
}
// ----------------------------------------------------------------------------

template <typename DATATYPE>
void write_field_elems(const Metadata& gmsh_options, const functionspace::CellColumns& function_space,
                       const Field& field, std::ostream& out) {
#if 1
    Log::debug() << "writing CellColumns field " << field.name() << "..." << std::endl;

    bool gather(gmsh_options.get<bool>("gather") && mpi::size() > 1);
    // unused: bool binary( !gmsh_options.get<bool>( "ascii" ) );
    idx_t nlev  = std::max<idx_t>(1, field.levels());
    idx_t ndata = std::min<idx_t>(function_space.nb_cells(), field.shape(0));
    idx_t nvars = std::max<idx_t>(1, field.variables());
    auto gidx   = array::make_view<gidx_t, 1>(function_space.cells().global_index());
    Field gidx_glb;
    Field field_glb;
    if (gather) {
        gidx_glb =
            function_space.createField<gidx_t>(option::name("gidx_glb") | option::levels(false) | option::global());
        function_space.gather(function_space.cells().global_index(), gidx_glb);
        gidx = array::make_view<gidx_t, 1>(gidx_glb);

        field_glb = function_space.createField(field, option::global());
        function_space.gather(field, field_glb);
        ndata = std::min<idx_t>(function_space.nb_cells_global(), field_glb.shape(0));
    }

    std::vector<int> lev = get_levels(nlev, gmsh_options);
    for (size_t ilev = 0; ilev < lev.size(); ++ilev) {
        int jlev = lev[ilev];
        if ((gather && mpi::rank() == 0) || !gather) {
            out << "$ElementData\n";
            out << "1\n";
            out << "\"" << field.name() << field_lev(field, jlev) << "\"\n";
            out << "1\n";
            out << field_time(field) << "\n";
            out << "4\n";
            out << field_step(field) << "\n";
            out << field_vars(nvars) << "\n";
            out << ndata << "\n";
            out << mpi::rank() << "\n";
            auto data = gather ? make_level_view<DATATYPE>(field_glb, ndata, jlev)
                               : make_level_view<DATATYPE>(field, ndata, jlev);
            write_level(out, gidx, data);
            out << "$EndElementData\n";
        }
    }
#endif
}

// ----------------------------------------------------------------------------
#if 0
template< typename DATA_TYPE >
void write_field_elems(const Metadata& gmsh_options, const FunctionSpace& function_space, const Field& field, std::ostream& out)
{
  Log::info() << "writing field " << field.name() << "..." << std::endl;
  bool gather( gmsh_options.get<bool>("gather") );
  bool binary( !gmsh_options.get<bool>("ascii") );
  int nlev = field.metadata().has("nb_levels") ? field.metadata().get<size_t>("nb_levels") : 1;
  int ndata = field.shape(0);
  int nvars = field.shape(1)/nlev;
  array::ArrayView<gidx_t,1    > gidx ( function_space.field( "glb_idx" ) );
  array::ArrayView<DATA_TYPE> data ( field );
  array::ArrayT<DATA_TYPE> field_glb_arr;
  array::ArrayT<gidx_t   > gidx_glb_arr;
  if( gather )
  {
    mpl::GatherScatter& fullgather = function_space.fullgather();
    ndata = fullgather.glb_dof();
    field_glb_arr.resize(ndata,field.shape(1));
    gidx_glb_arr.resize(ndata);
    array::ArrayView<DATA_TYPE> data_glb( field_glb_arr );
    array::ArrayView<gidx_t,1> gidx_glb( gidx_glb_arr );
    fullgather.gather( gidx, gidx_glb );
    fullgather.gather( data, data_glb );
    gidx = array::ArrayView<gidx_t,1>( gidx_glb_arr );
    data = data_glb;
  }

  double time = field.metadata().has("time") ? field.metadata().get<double>("time") : 0.;
  size_t step = field.metadata().has("step") ? field.metadata().get<size_t>("step") : 0 ;

  int nnodes = IndexView<int,2>( function_space.field("nodes") ).shape(1);

  for (int jlev=0; jlev<nlev; ++jlev)
  {
    char field_lev[6] = {0, 0, 0, 0, 0, 0};
    if( field.metadata().has("nb_levels") )
      std::sprintf(field_lev, "[%03d]",jlev);

    out << "$ElementNodeData\n";
    out << "1\n";
    out << "\"" << field.name() << field_lev << "\"\n";
    out << "1\n";
    out << time << "\n";
    out << "4\n";
    out << step << "\n";
    if     ( nvars == 1 ) out << nvars << "\n";
    else if( nvars <= 3 ) out << 3     << "\n";
    out << ndata << "\n";
    out << mpi::rank() << "\n";

    if( binary )
    {
      if( nvars == 1)
      {
        double value;
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out.write(reinterpret_cast<const char*>(&gidx(jelem)),sizeof(int));
          out.write(reinterpret_cast<const char*>(&nnodes),sizeof(int));
          for (size_t n=0; n<nnodes; ++n)
          {
            value = data(jelem,jlev);
            out.write(reinterpret_cast<const char*>(&value),sizeof(double));
          }
        }
      }
      else if( nvars <= 3 )
      {
        double value[3] = {0,0,0};
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
          {
            for( int v=0; v<nvars; ++v)
              value[v] = data(jelem,jlev*nvars+v);
            out.write(reinterpret_cast<const char*>(&value),sizeof(double)*3);
          }
        }
      }
      out <<"\n";
    }
    else
    {
      if( nvars == 1)
      {
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
            out << " " << data(jelem,jlev);
          out <<"\n";
        }
      }
      else if( nvars <= 3 )
      {
        std::vector<DATA_TYPE> data_vec(3,0.);
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
          {
            for( int v=0; v<nvars; ++v)
              data_vec[v] = data(jelem,jlev*nvars+v);
            for( int v=0; v<3; ++v)
              out << " " << data_vec[v];
          }
          out <<"\n";
        }
      }
    }
    out << "$EndElementNodeData\n";
  }
}
#endif

// ----------------------------------------------------------------------------

}  // end anonymous namespace

// ----------------------------------------------------------------------------
GmshIO::GmshIO() {
    // which field holds the Nodes
    options.set<std::string>("nodes", "xy");

    // Gather fields to one proc before writing
    options.set<bool>("gather", false);

    // Output of ghost nodes / elements
    options.set<bool>("ghost", false);

    // ASCII format (true) or binary (false)
    options.set<bool>("ascii", true);

    // Output of elements
    options.set<bool>("elements", true);

    // Output of edges
    options.set<bool>("edges", true);

    // Levels of fields to use
    options.set<std::vector<long>>("levels", std::vector<long>());
}

GmshIO::~GmshIO() = default;

Mesh GmshIO::read(const PathName& file_path) const {
    Mesh mesh;
    GmshIO::read(file_path, mesh);
    return mesh;
}

namespace {
mesh::ElementType* make_element_type(int type) {
    if (type == QUAD) {
        return new mesh::temporary::Quadrilateral();
    }
    if (type == TRIAG) {
        return new mesh::temporary::Triangle();
    }
    if (type == LINE) {
        return new mesh::temporary::Line();
    }
    throw_Exception("Element type not supported", Here());
}
}  // namespace

void GmshIO::read(const PathName& file_path, Mesh& mesh) const {
    std::ifstream file;
    file.open(file_path.localPath(), std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        throw_CantOpenFile(file_path);
    }

    std::string line;

    while (line != "$MeshFormat") {
        std::getline(file, line);
    }
    double version;
    int binary;
    int size_of_real;
    file >> version >> binary >> size_of_real;

    while (line != "$Nodes") {
        std::getline(file, line);
    }

    // Create nodes
    idx_t nb_nodes;
    file >> nb_nodes;

    mesh.nodes().resize(nb_nodes);

    mesh::Nodes& nodes = mesh.nodes();

    //nodes.add( Field( "xyz", array::make_datatype<double>(), array::make_shape( nb_nodes, 3 ) ) );

    //    array::ArrayView<double, 2> coords  = array::make_view<double, 2>( nodes.field( "xyz" ) );
    array::ArrayView<double, 2> xy      = array::make_view<double, 2>(nodes.xy());
    array::ArrayView<double, 2> lonlat  = array::make_view<double, 2>(nodes.lonlat());
    array::ArrayView<gidx_t, 1> glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());
    array::ArrayView<int, 1> part       = array::make_view<int, 1>(nodes.partition());
    array::ArrayView<int, 1> ghost      = array::make_view<int, 1>(nodes.ghost());

    std::map<int, int> glb_to_loc;
    int g;
    double x, y, z;
    double xyz[3];
    double xmax        = -std::numeric_limits<double>::max();
    double zmax        = -std::numeric_limits<double>::max();
    gidx_t max_glb_idx = 0;
    while (binary && file.peek() == '\n') {
        file.get();
    }
    for (idx_t n = 0; n < nb_nodes; ++n) {
        if (binary) {
            file.read(reinterpret_cast<char*>(&g), sizeof(int));
            file.read(reinterpret_cast<char*>(&xyz), sizeof(double) * 3);
            x = xyz[XX];
            y = xyz[YY];
            z = xyz[ZZ];
        }
        else {
            file >> g >> x >> y >> z;
        }
        glb_idx(n)     = g;
        xy(n, XX)      = x;
        xy(n, YY)      = y;
        lonlat(n, LON) = x;
        lonlat(n, LAT) = y;
        glb_to_loc[g]  = n;
        part(n)        = 0;
        ghost(n)       = 0;
        max_glb_idx    = std::max(max_glb_idx, static_cast<gidx_t>(g));
        xmax           = std::max(x, xmax);
        zmax           = std::max(z, zmax);
    }
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
    }

    int nb_elements = 0;

    while (line != "$Elements") {
        std::getline(file, line);
    }

    file >> nb_elements;

    if (binary) {
        while (file.peek() == '\n') {
            file.get();
        }
        int accounted_elems = 0;
        while (accounted_elems < nb_elements) {
            int header[3];
            int data[100];
            file.read(reinterpret_cast<char*>(&header), sizeof(int) * 3);

            int etype     = header[0];
            size_t netype = header[1];
            size_t ntags  = header[2];
            accounted_elems += netype;
            mesh::Elements* elements;
            if (etype == LINE) {
                size_t jtype = mesh.edges().add(make_element_type(etype), netype);
                elements     = &mesh.edges().elements(jtype);
            }
            else {
                size_t jtype = mesh.cells().add(make_element_type(etype), netype);
                elements     = &mesh.cells().elements(jtype);
            }

            size_t nnodes_per_elem            = elements->element_type().nb_nodes();
            mesh::BlockConnectivity& conn     = elements->node_connectivity();
            array::ArrayView<gidx_t, 1> egidx = array::make_view<gidx_t, 1>(elements->global_index());
            array::ArrayView<int, 1> epart    = array::make_view<int, 1>(elements->partition());

            size_t dsize = 1 + ntags + nnodes_per_elem;
            int part;
            for (size_t e = 0; e < netype; ++e) {
                file.read(reinterpret_cast<char*>(&data), sizeof(int) * dsize);
                part     = 0;
                egidx(e) = data[0];
                epart(e) = part;
                for (size_t n = 0; n < nnodes_per_elem; ++n) {
                    conn.set(e, n, glb_to_loc[data[1 + ntags + n]]);
                }
            }
        }
    }
    else {
        // Find out which element types are inside
        int position = file.tellg();
        std::vector<int> nb_etype(20, 0);
        int elements_max_glb_idx(0);
        int etype;
        for (int e = 0; e < nb_elements; ++e) {
            file >> g >> etype;
            std::getline(file, line);  // finish line
            ++nb_etype[etype];
            elements_max_glb_idx = std::max(elements_max_glb_idx, g);
        }

        // Allocate data structures for quads, triags, edges

        int nb_quads  = nb_etype[QUAD];
        int nb_triags = nb_etype[TRIAG];
        //int nb_edges  = nb_etype[LINE];

        mesh::Elements& quads  = mesh.cells().elements(mesh.cells().add(make_element_type(QUAD), nb_quads));
        mesh::Elements& triags = mesh.cells().elements(mesh.cells().add(make_element_type(TRIAG), nb_triags));
        //        mesh::Elements& edges  = mesh.edges().elements( mesh.edges().add( make_element_type( LINE ), nb_edges ) );

        mesh::BlockConnectivity& quad_nodes  = quads.node_connectivity();
        mesh::BlockConnectivity& triag_nodes = triags.node_connectivity();
        //        mesh::BlockConnectivity& edge_nodes  = edges.node_connectivity();

        array::ArrayView<gidx_t, 1> quad_glb_idx = array::make_view<gidx_t, 1>(quads.global_index());
        array::ArrayView<int, 1> quad_part       = array::make_view<int, 1>(quads.partition());

        array::ArrayView<gidx_t, 1> triag_glb_idx = array::make_view<gidx_t, 1>(triags.global_index());
        array::ArrayView<int, 1> triag_part       = array::make_view<int, 1>(triags.partition());

        //        array::ArrayView<gidx_t, 1> edge_glb_idx = array::make_view<gidx_t, 1>( edges.global_index() );
        //        array::ArrayView<int, 1> edge_part       = array::make_view<int, 1>( edges.partition() );

        // Now read all elements
        file.seekg(position, std::ios::beg);
        int gn0, gn1, gn2, gn3;
        int quad = 0, triag = 0;
        int ntags, tags[100];
        for (int e = 0; e < nb_elements; ++e) {
            file >> g >> etype >> ntags;
            for (int t = 0; t < ntags; ++t) {
                file >> tags[t];
            }
            int part = 0;
            if (ntags > 3) {
                part = std::max(part, *std::max_element(tags + 3, tags + ntags - 1));  // one positive, others negative
            }

            idx_t enodes[4] = {-1, -1, -1, -1};

            switch (etype) {
                case (QUAD):
                    file >> gn0 >> gn1 >> gn2 >> gn3;
                    quad_glb_idx(quad) = g;
                    quad_part(quad)    = part;
                    enodes[0]          = glb_to_loc[gn0];
                    enodes[1]          = glb_to_loc[gn1];
                    enodes[2]          = glb_to_loc[gn2];
                    enodes[3]          = glb_to_loc[gn3];
                    quad_nodes.set(quad, enodes);
                    ++quad;
                    break;
                case (TRIAG):
                    file >> gn0 >> gn1 >> gn2;
                    triag_glb_idx(triag) = g;
                    triag_part(triag)    = part;
                    enodes[0]            = glb_to_loc[gn0];
                    enodes[1]            = glb_to_loc[gn1];
                    enodes[2]            = glb_to_loc[gn2];
                    triag_nodes.set(triag, enodes);
                    ++triag;
                    break;
                case (LINE):
                    file >> gn0 >> gn1;
                    //                    edge_glb_idx( edge ) = g;
                    //                    edge_part( edge )    = part;
                    //                    enodes[0]            = glb_to_loc[gn0];
                    //                    enodes[1]            = glb_to_loc[gn1];
                    //                    edge_nodes.set( edge, enodes );
                    //                    ++edge;
                    break;
                case (POINT):
                    file >> gn0;
                    break;
                default:
                    std::cout << "etype " << etype << std::endl;
                    throw_Exception("ERROR: element type not supported", Here());
            }
        }
    }

    file.close();
}

void GmshIO::write(const Mesh& mesh, const PathName& file_path) const {
    int part           = mesh.metadata().has("part") ? mesh.metadata().get<size_t>("part") : mpi::rank();
    bool include_ghost = options.get<bool>("ghost") && options.get<bool>("elements");

    bool land_water_flag = options.has("water") || options.has("land");
    bool include_water   = options.getBool("water", false);
    bool include_land    = options.getBool("land", false);

    std::string nodes_field  = options.get<std::string>("nodes");
    const mesh::Nodes& nodes = mesh.nodes();

    const Field coords_field = nodes.field(nodes_field);
    array::ArrayT<double> dummy_double(1, 1);
    array::ArrayT<idx_t> dummy_idx(1, 1);
    bool coords_is_idx = coords_field.datatype().kind() == array::make_datatype<idx_t>().kind();
    auto coords        = array::make_view<const double, 2>(coords_is_idx ? dummy_double : coords_field.array());
    auto coords_idx    = array::make_view<const idx_t, 2>(coords_is_idx ? coords_field.array() : dummy_idx);

    auto glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());

    const idx_t surfdim = nodes.field(nodes_field).shape(1);  // nb of variables in coords

    bool include_patch = (surfdim == 3);


    ATLAS_ASSERT(surfdim == 2 || surfdim == 3);

    Log::debug() << "writing mesh to gmsh file " << file_path << std::endl;

    bool binary = !options.get<bool>("ascii");

    openmode mode = std::ios::out;
    if (binary) {
        mode = std::ios::out | std::ios::binary;
    }
    GmshFile file(file_path, mode, part);

    // Header
    if (binary) {
        write_header_binary(file);
    }
    else {
        write_header_ascii(file);
    }

    // Nodes
    const idx_t nb_nodes = nodes.size();
    file << "$Nodes\n";
    file << nb_nodes << "\n";
    double xyz[3] = {0., 0., 0.};
    for (idx_t n = 0; n < nb_nodes; ++n) {
        gidx_t g = glb_idx(n);

        if (coords_is_idx) {
            for (idx_t d = 0; d < surfdim; ++d) {
                xyz[d] = coords_idx(n, d);
            }
        }
        else {
            for (idx_t d = 0; d < surfdim; ++d) {
                xyz[d] = coords(n, d);
            }
        }

        if (binary) {
            file.write(reinterpret_cast<const char*>(&g), sizeof(gidx_t));
            file.write(reinterpret_cast<const char*>(&xyz), sizeof(double) * 3);
        }
        else {
            file << g << " " << xyz[XX] << " " << xyz[YY] << " " << xyz[ZZ] << "\n";
        }
    }
    if (binary) {
        file << "\n";
    }
    file << "$EndNodes\n";

    // Elements
    file << "$Elements\n";
    {
        std::vector<const mesh::HybridElements*> grouped_elements;
        if (options.get<bool>("elements")) {
            grouped_elements.push_back(&mesh.cells());
        }
        if (options.get<bool>("edges")) {
            grouped_elements.push_back(&mesh.edges());
        }

        idx_t nb_elements(0);
        for (const mesh::HybridElements* hybrid : grouped_elements) {
            nb_elements += hybrid->size();
            const auto hybrid_halo  = array::make_view<int, 1>(hybrid->halo());
            const auto hybrid_flags = array::make_view<int, 1>(hybrid->flags());

            auto include = [&](idx_t e) {
                auto topology = Topology::view(hybrid_flags(e));
                if (land_water_flag && not(include_water && include_land)) {
                    if (include_water && !topology.check(Topology::WATER)) {
                        return false;
                    }
                    if (include_land && !topology.check(Topology::LAND)) {
                        return false;
                    }
                }
                if (not include_ghost) {
                    if (topology.check(Topology::GHOST) || hybrid_halo(e)) {
                        return false;
                    }
                }
                if (not include_patch) {
                    if (topology.check(Topology::PATCH)) {
                        return false;
                    }
                }
                if (topology.check(Topology::INVALID)) {
                    return false;
                }
                return true;
            };

            for (idx_t e = 0; e < hybrid->size(); ++e) {
                if (not include(e)) {
                    --nb_elements;
                }
            }
        }

        file << nb_elements << "\n";

        for (const mesh::HybridElements* hybrid : grouped_elements) {
            for (idx_t etype = 0; etype < hybrid->nb_types(); ++etype) {
                const mesh::Elements& elements                   = hybrid->elements(etype);
                const mesh::ElementType& element_type            = elements.element_type();
                const mesh::BlockConnectivity& node_connectivity = elements.node_connectivity();
                size_t nb_nodes                                  = node_connectivity.cols();

                int gmsh_elem_type;
                if (element_type.name() == "Line") {
                    gmsh_elem_type = 1;
                }
                else if (element_type.name() == "Triangle") {
                    gmsh_elem_type = 2;
                }
                else if (element_type.name() == "Quadrilateral") {
                    gmsh_elem_type = 3;
                }
                else if (element_type.name() == "Pentagon") {
                    // Hack: treat as quadrilateral and ignore 5th point
                    gmsh_elem_type = 3;
                    nb_nodes       = 4;
                }
                else {
                    ATLAS_NOTIMPLEMENTED;
                }


                auto elems_glb_idx   = elements.view<gidx_t, 1>(elements.global_index());
                auto elems_partition = elements.view<int, 1>(elements.partition());

                auto elems_halo  = elements.view<int, 1>(elements.halo());
                auto elems_flags = elements.view<int, 1>(elements.flags());

                auto include = [&](idx_t e) {
                    auto topology = Topology::view(elems_flags(e));
                    if (land_water_flag && not(include_water && include_land)) {
                        if (include_water && !topology.check(Topology::WATER)) {
                            return false;
                        }
                        if (include_land && !topology.check(Topology::LAND)) {
                            return false;
                        }
                    }
                    if (not include_ghost) {
                        if (topology.check(Topology::GHOST) || elems_halo(e)) {
                            return false;
                        }
                    }
                    if (not include_patch) {
                        if (topology.check(Topology::PATCH)) {
                            return false;
                        }
                    }
                    if (topology.check(Topology::INVALID)) {
                        return false;
                    }
                    return true;
                };
                if (binary) {
                    idx_t nb_elems = elements.size();
                    if (!include_ghost) {
                        for (idx_t elem = 0; elem < elements.size(); ++elem) {
                            if (elems_halo(elem)) {
                                --nb_elems;
                            }
                        }
                    }

                    int header[3];
                    int data[9];
                    header[0] = gmsh_elem_type;
                    header[1] = nb_elems;
                    header[2] = 4;  // nb_tags
                    file.write(reinterpret_cast<const char*>(&header), sizeof(int) * 3);
                    data[1]         = 1;
                    data[2]         = 1;
                    data[3]         = 1;
                    size_t datasize = sizeof(int) * (5 + nb_nodes);
                    for (idx_t elem = 0; elem < nb_elems; ++elem) {
                        if (include_ghost || !elems_halo(elem)) {
                            data[0] = elems_glb_idx(elem);
                            data[4] = elems_partition(elem);
                            for (idx_t n = 0; n < nb_nodes; ++n) {
                                data[5 + n] = glb_idx(node_connectivity(elem, n));
                            }
                            file.write(reinterpret_cast<const char*>(&data), datasize);
                        }
                    }
                }
                else {
                    std::stringstream ss_elem_info;
                    ss_elem_info << " " << gmsh_elem_type << " 4 1 1 1 ";
                    std::string elem_info = ss_elem_info.str();
                    for (idx_t elem = 0; elem < elements.size(); ++elem) {
                        if (include(elem)) {
                            file << elems_glb_idx(elem) << elem_info << elems_partition(elem);
                            for (idx_t n = 0; n < nb_nodes; ++n) {
                                file << " " << glb_idx(node_connectivity(elem, n));
                            }
                            file << "\n";
                        }
                    }
                }
            }
        }
    }
    if (binary) {
        file << "\n";
    }
    file << "$EndElements\n";
    file << std::flush;

    // Optional mesh information file
    if (options.has("info") && options.get<bool>("info")) {
        PathName mesh_info(file_path);
        mesh_info = mesh_info.dirName() + "/" + mesh_info.baseName(false) + "_info.msh";

        //[next]  make NodesFunctionSpace accept const mesh
        functionspace::NodeColumns function_space(const_cast<Mesh&>(mesh));

        FieldSet fieldset;
        auto lat =
            array::make_view<double, 1>(fieldset.add(Field("lat", array::make_datatype<double>(), {nodes.size()})));
        auto lon =
            array::make_view<double, 1>(fieldset.add(Field("lon", array::make_datatype<double>(), {nodes.size()})));

        auto lonlat = array::make_view<double, 2>(nodes.lonlat());
        for (idx_t n = 0; n < nodes.size(); ++n) {
            lon(n) = lonlat(n, 0);
            lat(n) = lonlat(n, 1);
        }
        write(fieldset, function_space, mesh_info, std::ios_base::out);
        std::vector<std::string> extra_fields = {"partition", "water", "dual_volumes", "dual_delta_sph",
                                                 "ghost",     "halo",  "remote_index"};
        for (auto& f : extra_fields) {
            if (nodes.has_field(f)) {
                write(nodes.field(f), function_space, mesh_info, std::ios_base::app);
            }
        }


        //[next] if( mesh.has_function_space("edges") )
        //[next] {
        //[next]   FunctionSpace& edges = mesh.function_space( "edges" );

        //[next]   if (edges.has_field("dual_normals"))
        //[next]   {
        //[next] write(edges.field("dual_normals"),mesh_info,std::ios_base::app);
        //[next]   }

        //[next]   if (edges.has_field("skewness"))
        //[next]   {
        //[next]     write(edges.field("skewness"),mesh_info,std::ios_base::app);
        //[next]   }

        //[next]   if (edges.has_field("arc_length"))
        //[next]   {
        //[next]     write(edges.field("arc_length"),mesh_info,std::ios_base::app);
        //[next]   }
        //[next] }
    }
    file.close();
}

// ----------------------------------------------------------------------------
void GmshIO::write(const Field& field, const PathName& file_path, openmode mode) const {
    if (!field.functionspace()) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }

    else if (functionspace::NodeColumns(field.functionspace())) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }
    else if (functionspace::StructuredColumns(field.functionspace())) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }
    else if (functionspace::CellColumns(field.functionspace())) {
        FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, field.functionspace(), file_path, mode);
    }
    else {
        std::stringstream msg;
        msg << "Field [" << field.name() << "] has functionspace [" << field.functionspace().type()
            << "] but requires a [functionspace::NodeColumns "
            << "or functionspace::StructuredColumns]";

        throw_AssertionFailed(msg.str(), Here());
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const Field& field, const functionspace::NodeColumns& functionspace,
                            const PathName& file_path, openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}

// ----------------------------------------------------------------------------

void GmshIO::write_delegate(const Field& field, const functionspace::NoFunctionSpace& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}

// ----------------------------------------------------------------------------

void GmshIO::write_delegate(const Field& field, const functionspace::CellColumns& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const Field& field, const functionspace::StructuredColumns& functionspace,
                            const PathName& file_path, openmode mode) const {
    FieldSet fieldset;
    fieldset.add(field);
    write_delegate(fieldset, functionspace, file_path, mode);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::NodeColumns& functionspace,
                            const PathName& file_path, openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));
    if (binary) {
        mode |= std::ios_base::binary;
    }
    bool gather = options.has("gather") ? options.get<bool>("gather") : false;
    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    // Header
    if (is_new_file) {
        write_header_ascii(file);
    }

    // field::Fields
    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_nodes<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_nodes<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_nodes<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_nodes<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }
    file.close();
}

void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::NoFunctionSpace& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));
    if (binary) {
        mode |= std::ios_base::binary;
    }
    bool gather = options.has("gather") ? options.get<bool>("gather") : false;
    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    // Header
    if (is_new_file) {
        write_header_ascii(file);
    }

    // field::Fields
    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_nodes<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_nodes<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_nodes<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_nodes<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }
    file.close();
}

void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::CellColumns& functionspace,
                            const eckit::PathName& file_path, GmshIO::openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));
    if (binary) {
        mode |= std::ios_base::binary;
    }
    bool gather = options.has("gather") ? options.get<bool>("gather") : false;
    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    // Header
    if (is_new_file) {
        write_header_ascii(file);
    }

    // field::Fields
    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_elems<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_elems<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_elems<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_elems<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }
    file.close();
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write_delegate(const FieldSet& fieldset, const functionspace::StructuredColumns& functionspace,
                            const PathName& file_path, openmode mode) const {
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));

    if (binary) {
        mode |= std::ios_base::binary;
    }

    bool gather = options.has("gather") ? options.get<bool>("gather") : false;

    GmshFile file(file_path, mode, gather ? -1 : int(mpi::rank()));

    // Header
    if (is_new_file) {
        write_header_ascii(file);
    }

    // field::Fields
    for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
        const Field& field = fieldset[field_idx];
        Log::debug() << "writing field " << field.name() << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32()) {
            write_field_nodes<int>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64()) {
            write_field_nodes<long>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32()) {
            write_field_nodes<float>(options, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64()) {
            write_field_nodes<double>(options, functionspace, field, file);
        }

        file << std::flush;
    }

    file.close();
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write(const FieldSet& fieldset, const FunctionSpace& funcspace, const eckit::PathName& file_path,
                   openmode mode) const {
    if (functionspace::NodeColumns(funcspace)) {
        write_delegate(fieldset, functionspace::NodeColumns(funcspace), file_path, mode);
    }
    else if (functionspace::StructuredColumns(funcspace)) {
        write_delegate(fieldset, functionspace::StructuredColumns(funcspace), file_path, mode);
    }
    else if (functionspace::CellColumns(funcspace)) {
        write_delegate(fieldset, functionspace::CellColumns(funcspace), file_path, mode);
    }
    else if (not funcspace) {
        write_delegate(fieldset, functionspace::NoFunctionSpace(), file_path, mode);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void GmshIO::write(const Field& field, const FunctionSpace& funcspace, const eckit::PathName& file_path,
                   openmode mode) const {
    if (functionspace::NodeColumns(funcspace)) {
        write_delegate(field, functionspace::NodeColumns(funcspace), file_path, mode);
    }
    else if (functionspace::StructuredColumns(funcspace)) {
        write_delegate(field, functionspace::StructuredColumns(funcspace), file_path, mode);
    }
    else if (functionspace::CellColumns(funcspace)) {
        write_delegate(field, functionspace::CellColumns(funcspace), file_path, mode);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}
// ----------------------------------------------------------------------------

class GmshFortranInterface {
public:
    static Mesh::Implementation* atlas__Gmsh__read(GmshIO* This, char* file_path);
    static void atlas__Gmsh__write(GmshIO* This, Mesh::Implementation* mesh, char* file_path);
    static Mesh::Implementation* atlas__read_gmsh(char* file_path);
    static void atlas__write_gmsh_mesh(const Mesh::Implementation* mesh, char* file_path);
    static void atlas__write_gmsh_fieldset(const field::FieldSetImpl* fieldset,
                                           functionspace::FunctionSpaceImpl* functionspace, char* file_path, int mode);
    static void atlas__write_gmsh_field(const field::FieldImpl* field, functionspace::FunctionSpaceImpl* functionspace,
                                        char* file_path, int mode);
};

Mesh::Implementation* GmshFortranInterface::atlas__Gmsh__read(GmshIO* This, char* file_path) {
    Mesh::Implementation* m;
    {
        Mesh mesh = This->read(PathName(file_path));
        mesh.get()->attach();
        m = mesh.get();
    }
    m->detach();
    return m;
}

void GmshFortranInterface::atlas__Gmsh__write(GmshIO* This, Mesh::Implementation* mesh, char* file_path) {
    Mesh m(mesh);
    This->write(m, PathName(file_path));
}

Mesh::Implementation* GmshFortranInterface::atlas__read_gmsh(char* file_path) {
    Mesh::Implementation* m;
    {
        Mesh mesh = GmshIO().read(PathName(file_path));
        mesh.get()->attach();
        m = mesh.get();
    }
    m->detach();
    return m;
}

void GmshFortranInterface::atlas__write_gmsh_mesh(const Mesh::Implementation* mesh, char* file_path) {
    GmshIO writer;
    writer.write(mesh, PathName(file_path));
}

void GmshFortranInterface::atlas__write_gmsh_fieldset(const field::FieldSetImpl* fieldset,
                                                      functionspace::FunctionSpaceImpl* functionspace, char* file_path,
                                                      int /*mode*/) {
    GmshIO writer;
    writer.write(fieldset, functionspace, PathName(file_path));
}

void GmshFortranInterface::atlas__write_gmsh_field(const field::FieldImpl* field,
                                                   functionspace::FunctionSpaceImpl* functionspace, char* file_path,
                                                   int /*mode*/) {
    GmshIO writer;
    writer.write(field, functionspace, PathName(file_path));
}

extern "C" {

// ----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
// ----------------------------------------------------------------------------
GmshIO* atlas__Gmsh__new() {
    return new GmshIO();
}

void atlas__Gmsh__delete(GmshIO* This) {
    delete This;
}

Mesh::Implementation* atlas__Gmsh__read(GmshIO* This, char* file_path) {
    return GmshFortranInterface::atlas__Gmsh__read(This, file_path);
}

void atlas__Gmsh__write(GmshIO* This, Mesh::Implementation* mesh, char* file_path) {
    GmshFortranInterface::atlas__Gmsh__write(This, mesh, file_path);
}

Mesh::Implementation* atlas__read_gmsh(char* file_path) {
    return GmshFortranInterface::atlas__read_gmsh(file_path);
}

void atlas__write_gmsh_mesh(const Mesh::Implementation* mesh, char* file_path) {
    GmshFortranInterface::atlas__write_gmsh_mesh(mesh, file_path);
}

void atlas__write_gmsh_fieldset(const field::FieldSetImpl* fieldset, functionspace::FunctionSpaceImpl* functionspace,
                                char* file_path, int mode) {
    GmshFortranInterface::atlas__write_gmsh_fieldset(fieldset, functionspace, file_path, mode);
}

void atlas__write_gmsh_field(const field::FieldImpl* field, functionspace::FunctionSpaceImpl* functionspace,
                             char* file_path, int mode) {
    GmshFortranInterface::atlas__write_gmsh_field(field, functionspace, file_path, mode);
}
}
// ----------------------------------------------------------------------------

}  // namespace detail
}  // namespace output
}  // namespace atlas
