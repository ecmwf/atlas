
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "eckit/runtime/Context.h"
#include "eckit/config/Configurable.h"
#include "eckit/geometry/Point3.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/grid/partitioners/EqualRegionsPartitioner.h"
#include "atlas/grid/regular/Regular.h"
#include "atlas/grid/GridDistribution.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/generators/RegularMeshGenerator.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Parameters.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/runtime/Log.h"
#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/internals/Debug.h"
#include "atlas/runtime/Log.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

using namespace eckit;

using atlas::internals::Topology;

namespace atlas {
namespace mesh {
namespace generators {

RegularMeshGenerator::RegularMeshGenerator(const eckit::Parametrisation& p)
{
  configure_defaults();

  // options copied from Structured MeshGenerator
  size_t nb_parts;
  if( p.get("nb_parts",nb_parts) )
    options.set("nb_parts",nb_parts);

  size_t part;
  if( p.get("part",part) )
    options.set("part",part);

  std::string partitioner;
  if( p.get("partitioner",partitioner) )
  {
    if( not grid::partitioners::PartitionerFactory::has(partitioner) ) {
      Log::warning() << "Atlas does not have support for partitioner " << partitioner << ". "
                     << "Defaulting to use partitioner EqualRegions" << std::endl;
      partitioner = "EqualRegions";
    }
    options.set("partitioner",partitioner);
  }

  // options specifically for this MeshGenerator
  bool periodic_x;
  if( p.get("periodic_x",periodic_x) )
    options.set("periodic_x",periodic_x);

  bool periodic_y;
  if( p.get("periodic_y",periodic_y) )
    options.set("periodic_y",periodic_y);

  bool biperiodic;
  if( p.get("biperiodic",biperiodic) ) {
    options.set("periodic_x",biperiodic);
    options.set("periodic_y",biperiodic);
  }

  /***
  // original options from (global) structured meshgenerator, that don't seem to apply here:
  bool include_pole;
  if( p.get("include_pole",include_pole) )
    options.set("include_pole",include_pole);

  bool patch_pole;
  if( p.get("patch_pole",patch_pole) )
    options.set("patch_pole",patch_pole);

  bool unique_pole;
  if( p.get("unique_pole",unique_pole) )
    options.set("unique_pole",unique_pole);

  bool three_dimensional;
  if( p.get("three_dimensional",three_dimensional) || p.get("3d",three_dimensional) )
    options.set("3d",three_dimensional);

  double angle;
  if( p.get("angle",angle) )
    options.set("angle",angle);

  bool triangulate;
  if( p.get("triangulate",triangulate) )
    options.set("triangulate",triangulate);

  bool ghost_at_end;
  if( p.get("ghost_at_end",ghost_at_end) )
    options.set("ghost_at_end",ghost_at_end);

  ***/

}


void RegularMeshGenerator::configure_defaults()
{

  // This option sets number of parts the mesh will be split in
  options.set( "nb_parts", eckit::mpi::size() );

  // This option sets the part that will be generated
  options.set( "part", eckit::mpi::rank() );

  // This options sets the default partitioner
  std::string partitioner;
  if( grid::partitioners::PartitionerFactory::has("Trans") )
    partitioner = "Trans";
  else
    partitioner = "EqualRegions";
  options.set<std::string>("partitioner",partitioner);

  // Options for for periodic grids
  options.set<bool>("periodic_x",false);
  options.set<bool>("periodic_y",false);

  /***
  // original options from (global) structured meshgenerator, that don't seem to apply here:

  // This option creates a point at the pole when true
  options.set( "include_pole", false );

  // This option sets the part that will be generated
  options.set( "patch_pole", false );

  // This option disregards multiple poles in grid (e.g. lonlat up to poles) and connects elements
  // to the first node only. Note this option will only be looked at in case other option
  // "3d"==true
  options.set( "unique_pole", true );

  // This option creates elements that connect east to west at greenwich meridian
  // when true, instead of creating periodic ghost-points at east boundary when false
  options.set( "3d", false );

  // Experimental option. The result is a non-standard Reduced Gaussian Grid, with a ragged Greenwich line
  options.set("stagger", false );

  // This option sets the maximum angle deviation for a quadrilateral element
  // angle = 30  -->  minimises number of triangles
  // angle = 0   -->  maximises number of triangles
  options.set<double>("angle", 0. );

  options.set<bool>("triangulate", false );

  // This options moves the ghost points to the end
  options.set<bool>("ghost_at_end", true );

  ***/

}

void RegularMeshGenerator::generate(const grid::Grid& grid, Mesh& mesh ) const
{
    ASSERT(!mesh.generated());

  const grid::regular::Regular* rg = dynamic_cast<const grid::regular::Regular*>(&grid);
  if( !rg )
    throw eckit::BadCast("RegularMeshGenerator can only work with a Regular grid",Here());

  size_t nb_parts = options.get<size_t>("nb_parts");

  std::string partitioner_factory = "EqualRegions";
  options.get("partitioner",partitioner_factory);

  //if ( rg->nlat()%2 == 1 ) partitioner_factory = "EqualRegions"; // Odd number of latitudes
  //if ( nb_parts == 1 || eckit::mpi::size() == 1 ) partitioner_factory = "EqualRegions"; // Only one part --> Trans is slower

  grid::partitioners::Partitioner::Ptr partitioner( grid::partitioners::PartitionerFactory::build(partitioner_factory,grid,nb_parts) );
  grid::GridDistribution::Ptr distribution( partitioner->distribution() );
  generate( grid, *distribution, mesh );
}

void RegularMeshGenerator::hash(MD5& md5) const
{
    md5.add("RegularMeshGenerator");
#if ECKIT_MAJOR_VERSION >= 0 && ECKIT_MINOR_VERSION >= 13
    options.hash(md5);
#endif
}

void RegularMeshGenerator::generate(const grid::Grid& grid, const grid::GridDistribution& distribution, Mesh& mesh ) const
{
  const grid::regular::Regular* rg = dynamic_cast<const grid::regular::Regular*>(&grid);
  if( !rg )
    throw eckit::BadCast("Grid could not be cast to a Regular",Here());

  ASSERT(!mesh.generated());

  if( grid.npts() != distribution.partition().size() )
  {
    std::stringstream msg;
    msg << "Number of points in grid ("<<grid.npts()<<") different from "
           "number of points in grid distribution ("<<distribution.partition().size()<<")";
    throw eckit::AssertionFailed(msg.str(),Here());
  }

  int mypart   = options.get<size_t>("part");

  // clone some grid properties
  mesh.setProjection(rg->projection());

  generate_mesh(*rg,distribution,mesh);
}

void RegularMeshGenerator::generate_mesh(
    const grid::regular::Regular& rg,
    const std::vector<int>& parts,
    //const Region& region,
    Mesh& mesh ) const
{
  int mypart = options.get<size_t>("part");
  int nparts = options.get<size_t>("nb_parts");
  int nx=rg.nlonmin();
  int ny=rg.nlat();

  bool periodic_x = options.get<bool>("periodic_x");
  bool periodic_y = options.get<bool>("periodic_y");

  // for asynchronous output
#if DEBUG_OUTPUT
  sleep( mypart );
#endif

  // this function should do the following:
  // - define nodes with
  //      mesh.nodes().resize(nnodes);
  //      mesh::Nodes& nodes = mesh.nodes();
  //    following properties should be defined:
  //      array::ArrayView<double,2> lonlat        ( nodes.lonlat() );
  //      array::ArrayView<gidx_t,1> glb_idx       ( nodes.global_index() );
  //      array::ArrayView<int,   1> part          ( nodes.partition() );
  //      array::ArrayView<int,   1> ghost         ( nodes.ghost() );
  //      array::ArrayView<int,   1> flags         ( nodes.field("flags") );
  // - define cells (only quadrilaterals for now) with
  //      mesh.cells().add( new mesh::temporary::Quadrilateral(), nquads  );
  //    further define cells with
  //      array::ArrayView<gidx_t,1> cells_glb_idx( mesh.cells().global_index() );
  //      array::ArrayView<int,1>    cells_part(    mesh.cells().partition() );
  // - define connectivity with
  //      mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();
  //      node_connectivity.set( jcell, quad_nodes );
  //    where quad_nodes is a 4-element integer array containing the LOCAL indices of the nodes

  // Start with calculating number of quadrilaterals
  // The rule do determine if a cell belongs to a proc is the following: if the lowerleft corner of the cell belongs to that proc.
  // so we loop over all gridpoints, select those that belong to the proc, and determine the number of cells
  int ii_glb;  // global index
  int ncells;

  // vector of local indices: necessary for remote indices of ghost nodes
  std::vector<int> local_idx(rg.npts(),-1);
  std::vector<int> current_idx(nparts,0);  // index counter for each proc

  // determine rectangle (ix_min:ix_max) x (iy_min:iy_max) surrounding the nodes on this processor
  int ix_min, ix_max, iy_min, iy_max, ix_glb, iy_glb, ix, iy;
  int nnodes_nonghost, nnodes, nnodes_mx;    // number of nodes: non-ghost; total;  inside surrounding rectangle
  int ixr, iyu;  // indices of point to the right and above
  int nnodes_SR, ii;

  // loop over all points to determine local indices and surroundig rectangle
  ix_min=nx+1;ix_max=0;iy_min=ny+1;iy_max=0;
  nnodes_nonghost=0;

  ii_glb=0;
  for (iy=0;iy<ny; iy++) {
    for (ix=0;ix<nx;ix++) {
      local_idx[ii_glb]=current_idx[parts[ii_glb]]++;    // store local index on the local proc of this point
      if ( parts[ii_glb] == mypart )
      {
        ++nnodes_nonghost;  // non-ghost node: belongs to this part
        ix_min=std::min(ix_min,ix);
        ix_max=std::max(ix_max,ix);
        iy_min=std::min(iy_min,iy);
        iy_max=std::max(iy_max,iy);
      }
      ++ii_glb;                                         // global index
    }
  }

  // add one row/column for ghost nodes (which include periodicity points)
  ix_max=ix_max+1;
  iy_max=iy_max+1;

#if DEBUG_OUTPUT_DETAIL
        std::cout << "[" << mypart << "] : " << "SR = " << ix_min << ":" << ix_max << " x " << iy_min << ":" << iy_max << std::endl;
#endif

  // dimensions of surrounding rectangle (SR)
  int nxl=ix_max-ix_min+1;
  int nyl=iy_max-iy_min+1;

  // upper estimate for number of nodes
  nnodes_SR=nxl*nyl;

  // partitions and local indices in SR
  std::vector<int> parts_SR(nnodes_SR,-1);
  std::vector<int> local_idx_SR(nnodes_SR,-1);
  std::vector<bool> is_ghost_SR(nnodes_SR,true);
  ii=0;                                    // index inside SR
  for (iy=0; iy<nyl; iy++) {
    iy_glb=(iy_min+iy);                    // global y-index
    for (ix=0; ix<nxl; ix++) {
      ix_glb=(ix_min+ix);                  // global x-index
      is_ghost_SR[ii]=! ((parts_SR[ii]==mypart) && ix<nxl-1 && iy<nyl-1);
      if (ix_glb<nx && iy_glb<ny ) {
        ii_glb=(iy_glb)*nx+ix_glb;    // global index
        parts_SR[ii]=parts[ii_glb];
        local_idx_SR[ii]=local_idx[ii_glb];
        is_ghost_SR[ii]=! ((parts_SR[ii]==mypart) && ix<nxl-1 && iy<nyl-1);
      } else if (ix_glb==nx && iy_glb<ny) {
        // take properties from the point to the left
        parts_SR[ii]=parts[iy_glb*nx+ix_glb-1];
        local_idx_SR[ii]=-1;
        is_ghost_SR[ii]=true;
      } else if (iy_glb==ny && ix_glb<nx ) {
        // take properties from the point below
        parts_SR[ii]=parts[(iy_glb-1)*nx+ix_glb];
        local_idx_SR[ii]=-1;
        is_ghost_SR[ii]=true;
      } else {
        // take properties from the point belowleft
        parts_SR[ii]=parts[(iy_glb-1)*nx+ix_glb-1];
        local_idx_SR[ii]=-1;
        is_ghost_SR[ii]=true;
      }
      ++ii;
    }
  }

#if DEBUG_OUTPUT_DETAIL
  std::cout << "[" << mypart << "] : " << "parts_SR = "; for (ii=0;ii<nnodes_SR;ii++) std::cout << parts_SR[ii] << ","; std::cout << std::endl;
  std::cout << "[" << mypart << "] : " << "local_idx_SR = "; for (ii=0;ii<nnodes_SR;ii++) std::cout << local_idx_SR[ii] << ","; std::cout << std::endl;
  std::cout << "[" << mypart << "] : " << "is_ghost_SR = "; for (ii=0;ii<nnodes_SR;ii++) std::cout << is_ghost_SR[ii] << ","; std::cout << std::endl;
#endif

  // vectors marking nodes that are necessary for this proc's cells
  std::vector<bool> is_node_SR(nnodes_SR,false);

  // determine number of cells and number of nodes
  nnodes=0;
  ncells=0;
  for (iy=0; iy<nyl-1; iy++) {      // don't loop into ghost/periodicity row
    for (ix=0; ix<nxl-1; ix++) {    // don't loop into ghost/periodicity column
      ii=iy*nxl+ix;
      if ( ! is_ghost_SR[ii] ) {
        // mark this node as being used
        if (! is_node_SR[ii]) {
          ++nnodes;
          is_node_SR[ii]=true;
        }
        // check if this node is the lowerleft corner of a new cell
        if ( (ix_min+ix<nx-1 || periodic_x) && (iy_min+iy<ny-1 || periodic_y) ) {
          ++ncells;
          // mark lowerright corner
          ii=iy*nxl+ix+1;
          if (! is_node_SR[ii]) {
            ++nnodes;
            is_node_SR[ii]=true;
          }
          // mark upperleft corner
          ii=(iy+1)*nxl+ix;
          if (! is_node_SR[ii]) {
            ++nnodes;
            is_node_SR[ii]=true;
          }
          // mark upperright corner
          ii=(iy+1)*nxl+ix+1;
          if (! is_node_SR[ii]) {
            ++nnodes;
            is_node_SR[ii]=true;
          }
        }
        // periodic points are always needed, even if they don't belong to a cell
        ii=iy*nxl+ix+1;
        if ( periodic_x && ix_min+ix==nx-1 && ! is_node_SR[ii] ) {
          ++nnodes;
          is_node_SR[ii]=true;
        }
        ii=(iy+1)*nxl+ix;
        if ( periodic_y && iy_min+iy==ny-1 && ! is_node_SR[ii]  ) {
          ++nnodes;
          is_node_SR[ii]=true;
        }
        ii=(iy+1)*nxl+ix+1;
        if ( periodic_x && periodic_y && ix_min+ix==nx-1 && iy_min+iy==ny-1 && ! is_node_SR[ii]  ) {
          ++nnodes;
          is_node_SR[ii]=true;
        }

      }
    }
  }

#if DEBUG_OUTPUT_DETAIL
  std::cout << "[" << mypart << "] : " << "nnodes = " << nnodes << std::endl;
  std::cout << "[" << mypart << "] : " << "is_node_SR = "; for (ii=0;ii<nnodes_SR;ii++) std::cout << is_node_SR[ii] << ","; std::cout << std::endl;
#endif

  // define nodes and associated properties
  mesh.nodes().resize(nnodes);
  mesh::Nodes& nodes = mesh.nodes();
  array::ArrayView<double,2> lonlat        ( nodes.lonlat() );
  array::ArrayView<double,2> geolonlat     ( nodes.geolonlat() );
  array::ArrayView<gidx_t,1> glb_idx       ( nodes.global_index() );
  array::ArrayView<int,   1> remote_idx    ( nodes.remote_index() );
  array::ArrayView<int,   1> part          ( nodes.partition() );
  array::ArrayView<int,   1> ghost         ( nodes.ghost() );
  array::ArrayView<int,   1> flags         ( nodes.field("flags") );

  // define cells and associated properties
  mesh.cells().add( new mesh::temporary::Quadrilateral(), ncells );
  int quad_begin  = mesh.cells().elements(0).begin();
  array::ArrayView<int,1>    cells_part(    mesh.cells().partition() );
  mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();

  int quad_nodes[4];
  int jcell=quad_begin;
  int inode, inode_nonghost, inode_ghost;
  eckit::geometry::LLPoint2 llg;

  // global indices for periodicity points
  inode=nx*ny;
  std::vector<int> glb_idx_px(ny+1,-1);
  std::vector<int> glb_idx_py(nx+1,-1);
  if (periodic_x) {
    for (iy=0;iy<ny+(periodic_y?1:0);iy++) {
      glb_idx_px[iy]=inode++;
    }
  }
  if (periodic_y) {
    for (ix=0;ix<nx;ix++) {
      glb_idx_py[ix]=inode++;
    }
  }

  // loop over nodes and set properties
  ii=0;
  inode_nonghost=0;
  inode_ghost=nnodes_nonghost;  // ghost nodes start counting after nonghost nodes
  for (iy=0; iy<nyl; iy++) {
    for (ix=0; ix<nxl; ix++) {
      // node properties
      if ( is_node_SR[ii] ) {
        // set node counter
        if ( is_ghost_SR[ii] ) {
          inode=inode_ghost++;
        } else {
          inode=inode_nonghost++;
        }
        // global index
        ix_glb=(ix_min+ix);      // don't take modulus here: periodicity points have their own global index.
        iy_glb=(iy_min+iy);
        if ( ix_glb < nx && iy_glb < ny ) {
          ii_glb=iy_glb*nx+ix_glb;  // no periodic point
        } else {
          if ( ix_glb==nx ) {
            // periodicity point in x-direction
            ii_glb=glb_idx_px[iy_glb];
          } else {

            // periodicity point in x-direction
            ii_glb=glb_idx_py[ix_glb];
          }
        }
        glb_idx(inode) = ii_glb+1;  // starting from 1
        // grid coordinates
        double xy[2];
        if (iy_glb<ny) {
          // normal calculation
          rg.lonlat((size_t) iy_glb,(size_t) ix_glb,xy);
        } else {
          // for periodic_y grids, iy_glb==ny lies outside the range of latitudes in the Structured grid...
          // so we extrapolate from two other points -- this is okay for regular grids with uniform spacing.
          double xy1[2], xy2[2];
          rg.lonlat((size_t) iy_glb-1,(size_t) ix_glb,xy1);
          rg.lonlat((size_t) iy_glb-2,(size_t) ix_glb,xy2);
          xy[0]=2*xy1[0]-xy2[0];
          xy[1]=2*xy1[1]-xy2[1];
        }
        lonlat(inode,internals::LON) = xy[internals::LON];
        lonlat(inode,internals::LAT) = xy[internals::LAT];

        // geographic coordinates by using projection
        llg=rg.projection()->coords2lonlat(eckit::geometry::Point2(xy[0],xy[1]));
        geolonlat(inode,internals::LON) = llg[internals::LON];
        geolonlat(inode,internals::LAT) = llg[internals::LAT];

        // part
        part(inode) = parts_SR[ii];
        // ghost nodes
        ghost(inode)=is_ghost_SR[ii];
        // flags
        Topology::reset(flags(inode));
        if ( ghost(inode) ) {
          Topology::set(flags(inode),Topology::GHOST);
          remote_idx(inode)=local_idx_SR[ii];
          // change local index -- required for cells
          local_idx_SR[ii]=inode;
        } else {
          remote_idx(inode)=-1;
        }

#if DEBUG_OUTPUT_DETAIL
        std::cout << "[" << mypart << "] : " << "New node " << "\n\t";
        std::cout << "[" << mypart << "] : "  << "\tinode=" << inode << "; ix_glb=" << ix_glb << "; iy_glb=" << iy_glb << "; glb_idx=" << ii_glb << std::endl;
        std::cout << "[" << mypart << "] : "  << "\tglon=" << geolonlat(inode,0) << "; glat=" << geolonlat(inode,1) << "; glb_idx=" << glb_idx(inode) << std::endl;
#endif
      }
      ++ii;
    }
  }


  // loop over nodes and define cells
  for (iy=0; iy<nyl-1; iy++) {      // don't loop into ghost/periodicity row
    for (ix=0; ix<nxl-1; ix++) {    // don't loop into ghost/periodicity column
      ii=iy*nxl+ix;
      if ( ! is_ghost_SR[ii] ) {
        if ( (ix_min+ix<nx-1 || periodic_x) && (iy_min+iy<ny-1 || periodic_y) ) {
          // define cell corners (local indices)
          quad_nodes[0]=local_idx_SR[ii];
          quad_nodes[1]=local_idx_SR[iy*nxl+ix+1];      // point to the right
          quad_nodes[2]=local_idx_SR[(iy+1)*nxl+ix+1];  // point above right
          quad_nodes[3]=local_idx_SR[(iy+1)*nxl+ix];    // point above
          node_connectivity.set( jcell, quad_nodes );
          cells_part(jcell)    = mypart;
#if DEBUG_OUTPUT_DETAIL
          std::cout << "[" << mypart << "] : "  << "New quad " << jcell << "\n\t";
          std::cout << "[" << mypart << "] : "  << quad_nodes[0] << "," << quad_nodes[1] << "," << quad_nodes[2] << "," << quad_nodes[3] << std::endl;
#endif
          ++jcell;
        }
      }
    }
  }

#if DEBUG_OUTPUT
  // list nodes
  for (inode=0;inode<nnodes;inode++) {
    std::cout << "[" << mypart << "] : " << " node " << inode << ": ghost = " << ghost(inode) << ", glb_idx = " << glb_idx(inode)-1
      << ", part = " << part(inode) << ", lon = " << lonlat(inode,0) << ", lat = " << lonlat(inode,1)
      << ", remote_idx = " << remote_idx(inode)
      << std::endl;
  }

  int * cell_nodes;
  for (jcell=0;jcell<ncells;jcell++) {
    std::cout << "[" << mypart << "] : " << " cell " << jcell << ": " << node_connectivity(jcell,0) << ","
      << node_connectivity(jcell,1) << ","<< node_connectivity(jcell,2) << ","<< node_connectivity(jcell,3) << std::endl;
  }
#endif

  generate_global_element_numbering( mesh );

  nodes.metadata().set("parallel",true);

}


void RegularMeshGenerator::generate_global_element_numbering( Mesh& mesh ) const
{
  int loc_nb_elems = mesh.cells().size();
  std::vector<int> elem_counts( eckit::mpi::size() );
  std::vector<int> elem_displs( eckit::mpi::size() );

  ECKIT_MPI_CHECK_RESULT(
        MPI_Allgather( &loc_nb_elems, 1, MPI_INT,
                       elem_counts.data(), 1, MPI_INT, eckit::mpi::comm()) );
  elem_displs.at(0) = 0;
  for(size_t jpart = 1; jpart < eckit::mpi::size(); ++jpart)
  {
    elem_displs.at(jpart) = elem_displs.at(jpart-1) + elem_counts.at(jpart-1);
  }

  gidx_t gid = 1+elem_displs.at( eckit::mpi::rank() );

  array::ArrayView<gidx_t,1> glb_idx( mesh.cells().global_index() );

  for( size_t jelem=0; jelem<mesh.cells().size(); ++jelem )
  {
    glb_idx(jelem) = gid++;
  }
}

namespace {
static MeshGeneratorBuilder< RegularMeshGenerator > __RegularMeshGenerator("RegularMeshGenerator");
}

} // namespace generators
} // namespace mesh
} // namespace atlas
