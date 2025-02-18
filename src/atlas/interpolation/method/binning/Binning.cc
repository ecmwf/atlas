/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/config/LocalConfiguration.h"
#include "eckit/linalg/Triplet.h"
#include "eckit/linalg/SparseMatrix.h"
#include "eckit/mpi/Comm.h"

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/binning/Binning.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/mesh.h"
#include "atlas/mesh/actions/GetCubedSphereNodalArea.h"
#include "atlas/runtime/Trace.h"



namespace atlas {
namespace interpolation {
namespace method {


namespace {

MethodBuilder<Binning> __builder("binning");

}


Binning::Binning(const Config& config) : Method(config) {
  const auto* conf = dynamic_cast<const eckit::LocalConfiguration*>(&config);
  ATLAS_ASSERT(conf, "config must be derived from eckit::LocalConfiguration");
  interpAncillaryScheme_ = conf->getSubConfiguration("scheme");
  // enabling or disabling the adjoint operation
  adjoint_ = conf->getBool("adjoint", false);
  // enabling or disabling the halo exchange
  allow_halo_exchange_ = conf->getBool("halo_exchange", true);
}


void Binning::do_setup(const Grid& source,
                       const Grid& target,
                       const Cache&) {
  ATLAS_NOTIMPLEMENTED;
}

void Binning::do_setup(const FunctionSpace& source,
                       const FunctionSpace& target,
                       const Cache&) {
  ATLAS_NOTIMPLEMENTED;
}


void Binning::do_setup(const FunctionSpace& source,
                       const FunctionSpace& target) {
  ATLAS_TRACE("atlas::interpolation::method::Binning::do_setup()");

  using Scalar  = eckit::linalg::Scalar;
  using Index   = eckit::linalg::Index;
  using Triplet = eckit::linalg::Triplet;
  using SMatrix = eckit::linalg::SparseMatrix;

  source_ = source;
  target_ = target;

  if (target_.size() == 0) {
    return;
  }

  // note that the 'source' grid for the low-to-high regridding (interpolation)
  // is the 'target' grid for high-to-low regridding (binning) and
  // the 'target' grid for the low-to-high regridding (interpolation) is the
  // 'source' grid for the for high-to-low regridding (binning)
  const auto& fs_source_interp = target_;
  const auto& fs_target_interp = source_;

  const auto interp = Interpolation(
    interpAncillaryScheme_, fs_source_interp, fs_target_interp);
  auto smx_interp_cache = interpolation::MatrixCache(interp);

  auto smx_interp = smx_interp_cache.matrix();

  auto eckit_smx_interp = make_non_owning_eckit_sparse_matrix(smx_interp);
  SMatrix smx_interp_tr(eckit_smx_interp); // copy
  smx_interp_tr.transpose(); // transpose the copy in-place

  const auto rows_tamx = smx_interp_tr.rows();
  const auto cols_tamx = smx_interp_tr.cols();

  const Scalar* ptr_tamx_data = smx_interp_tr.data();
  const Index*  ptr_tamx_idxs_col = smx_interp_tr.inner();
  const Index*  ptr_tamx_o = smx_interp_tr.outer();

  // diagonal of 'area weights matrix', W
  auto ds_aweights = getAreaWeights(source_);

  auto smx_binning_els = std::vector<Triplet>{};
  size_t idx_row_next = 0;

  for (size_t idx_row = 0; idx_row < rows_tamx; ++idx_row) {
    idx_row_next = (idx_row+1);
    // start of the indexes associated with the row 'i'
    size_t lbound = ptr_tamx_o[idx_row];
    // start of the indexes associated with the row 'i+1'
    size_t ubound = ptr_tamx_o[idx_row_next];

    if (lbound == ubound) {
      continue;
    }

    double sum_row = 0;
    for (size_t i = lbound; i < ubound; ++i) {
      sum_row += (ptr_tamx_data[i] * ds_aweights.at(ptr_tamx_idxs_col[i]));
    }

    // normalization factor
    double nfactor = 1/sum_row;

    for (size_t i = lbound; i < ubound; ++i) {
      // evaluating the non-zero elements of the binning matrix
      smx_binning_els.emplace_back(
        idx_row, ptr_tamx_idxs_col[i],
        (nfactor * (ptr_tamx_data[i] * ds_aweights.at(ptr_tamx_idxs_col[i]))));
    }
  }

  // 'binning matrix' (sparse matrix), B = N A^T W
  setMatrix(rows_tamx, cols_tamx, smx_binning_els);
}


void Binning::print(std::ostream&) const {
  ATLAS_NOTIMPLEMENTED;
}


std::vector<double> Binning::getAreaWeights(const FunctionSpace& fspace) const {
  // diagonal of 'area weights matrix', W
  std::vector<double> ds_aweights;

  bool is_cubed_sphere {false};
  if (auto csfs = functionspace::NodeColumns(fspace)) {
    if (CubedSphereGrid(csfs.mesh().grid())) {
      is_cubed_sphere = true;
    }
  }

  if (is_cubed_sphere) { 

    const auto csfs = functionspace::NodeColumns(fspace); 
    auto csmesh = csfs.mesh();

    // areas of the cells (geographic coord. system)
    auto gcell_areas = mesh::actions::GetCubedSphereNodalArea()(csmesh);
    auto gcell_areas_view = array::make_view<double, 1>(gcell_areas);

    auto is_ghost = array::make_view<int, 1>(csfs.ghost());

    double total_area {0.};
    for (idx_t i = 0; i < csfs.size(); i++) {
      if (!is_ghost[i]) {
        total_area += gcell_areas_view(i);
      }
    }
    eckit::mpi::comm().allReduceInPlace(total_area, eckit::mpi::Operation::SUM);

    double aweight_temp {0.};
    for (idx_t i = 0; i < csfs.size(); i++) {
      if (!is_ghost[i]) {
        aweight_temp = gcell_areas_view(i)/total_area;
        ds_aweights.emplace_back(aweight_temp);
      }
    }
  
  } else {

    // area weights (default)
    ds_aweights.assign(fspace.size(), 1.);

  }

  return ds_aweights;
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
