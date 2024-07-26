/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/structured/RegionalLinear2D.h"

#include <limits>

#include "atlas/array.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace {
MethodBuilder<RegionalLinear2D> __builder("regional-linear-2d");
}

void RegionalLinear2D::print(std::ostream&) const { ATLAS_NOTIMPLEMENTED; }

void RegionalLinear2D::do_setup(const Grid& source, const Grid& target,
                                const Cache&) {
  ATLAS_NOTIMPLEMENTED;
}

void RegionalLinear2D::do_setup(const FunctionSpace& source,
                                const FunctionSpace& target) {
  ATLAS_TRACE("interpolation::method::RegionalLinear2D::do_setup");
  source_ = source;
  target_ = target;

  if (target_.size() == 0) {
    return;
  }
  ASSERT(source_.type() == "StructuredColumns");

  // Get grid parameters
  const functionspace::StructuredColumns sourceFs(source_);
  const RegularGrid sourceGrid(sourceFs.grid());
  const Projection & sourceProj = sourceGrid.projection();
  const size_t sourceNx = sourceGrid.nx();
  const size_t sourceNy = sourceGrid.ny();
  const double sourceDx = sourceGrid.dx();
  const double sourceDy = std::abs(sourceGrid.y(1)-sourceGrid.y(0));
  const bool reversedY = sourceGrid.y(1) < sourceGrid.y(0);

  // Check grid regularity in y direction
  for (size_t sourceJ = 0; sourceJ < sourceNy-1; ++sourceJ) {
    if (reversedY) {
      ASSERT(std::abs(sourceGrid.y(sourceJ)-sourceGrid.y(sourceJ+1)-sourceDy) < 1.0e-12*sourceDy);
    } else {
      ASSERT(std::abs(sourceGrid.y(sourceJ+1)-sourceGrid.y(sourceJ)-sourceDy) < 1.0e-12*sourceDy);
    }
  }

  // Source grid indices
  const Field sourceFieldIndexI = sourceFs.index_i();
  const Field sourceFieldIndexJ = sourceFs.index_j();
  const auto sourceIndexIView = array::make_view<int, 1>(sourceFieldIndexI);
  const auto sourceIndexJView = array::make_view<int, 1>(sourceFieldIndexJ);
  sourceSize_ = sourceFs.size();

  // Destination grid size
  targetSize_ = target_.size();

  // Ghost points
  const auto sourceGhostView = array::make_view<int, 1>(sourceFs.ghost());
  const auto targetGhostView = array::make_view<int, 1>(target_.ghost());

  // Define reduced grid horizontal distribution
  std::vector<int> mpiTask(sourceNx*sourceNy, 0);
  for (size_t sourceJnode = 0; sourceJnode < sourceSize_; ++sourceJnode) {
    if (sourceGhostView(sourceJnode) == 0) {
      mpiTask[(sourceIndexIView(sourceJnode)-1)*sourceNy+sourceIndexJView(sourceJnode)-1] = comm_.rank();
    }
  }
  comm_.allReduceInPlace(mpiTask.begin(), mpiTask.end(), eckit::mpi::sum());

  // Define local tree on destination grid
  std::vector<Point3> targetPoints;
  std::vector<size_t> targetIndices;
  const auto targetLonLatView = array::make_view<double, 2>(target_.lonlat());
  for (size_t targetJnode = 0; targetJnode < targetSize_; ++targetJnode) {
    PointLonLat p({targetLonLatView(targetJnode, 0), targetLonLatView(targetJnode, 1)});
    sourceProj.lonlat2xy(p);
    targetPoints.push_back(Point3(p[0], p[1], 0.0));
    targetIndices.push_back(targetJnode);
  }
  util::IndexKDTree targetTree;
  if (targetSize_ > 0) {
    targetTree.build(targetPoints, targetIndices);
  }
  const double radius = std::sqrt(sourceDx*sourceDx+sourceDy*sourceDy);

  // Delta for colocation
  const double eps = 1.0e-8;

  // RecvCounts and received points list
  targetRecvCounts_.resize(comm_.size());
  std::fill(targetRecvCounts_.begin(), targetRecvCounts_.end(), 0);
  std::vector<int> targetRecvPointsList;
  for (size_t sourceJ = 0; sourceJ < sourceNy; ++sourceJ) {
    double yMin, yMax;
    if (reversedY) {
      yMin = sourceJ < sourceNy-1 ? sourceGrid.y(sourceJ+1)-eps : -std::numeric_limits<double>::max();
      yMax = sourceJ > 0 ? sourceGrid.y(sourceJ-1)+eps : std::numeric_limits<double>::max();
    } else {
      yMin = sourceJ > 0 ? sourceGrid.y(sourceJ-1)-eps : -std::numeric_limits<double>::max();
      yMax = sourceJ < sourceNy-1 ? sourceGrid.y(sourceJ+1)+eps : std::numeric_limits<double>::max();
    }
    for (size_t sourceI = 0; sourceI < sourceNx; ++sourceI) {
      const double xMin = sourceI > 0 ? sourceGrid.x(sourceI-1)-eps : -std::numeric_limits<double>::max();
      const double xMax = sourceI < sourceNx-1 ? sourceGrid.x(sourceI+1)+eps :
        std::numeric_limits<double>::max();

      bool pointsNeeded = false;
      if (targetSize_ > 0) {
        const Point3 p(sourceGrid.x(sourceI), sourceGrid.y(sourceJ), 0.0);
        const auto list = targetTree.closestPointsWithinRadius(p, radius);
        for (const auto & item : list) {
          const PointXYZ targetPoint = item.point();
          const size_t targetJnode = item.payload();
          if (targetGhostView(targetJnode) == 0) {
            const bool inX = (xMin <= targetPoint[0] && targetPoint[0] <= xMax);
            const bool inY = (yMin <= targetPoint[1] && targetPoint[1] <= yMax);
            if (inX && inY) {
              pointsNeeded = true;
              break;
            }
          }
        }
      }
      if (pointsNeeded) {
        ++targetRecvCounts_[mpiTask[sourceI*sourceNy+sourceJ]];
        targetRecvPointsList.push_back(sourceI*sourceNy+sourceJ);
      }
    }
  }

  // Buffer size
  targetRecvSize_ = targetRecvPointsList.size();

  if (targetRecvSize_ > 0) {
    // RecvDispls
    targetRecvDispls_.push_back(0);
    for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
      targetRecvDispls_.push_back(targetRecvDispls_[jt]+targetRecvCounts_[jt]);
    }

    // Allgather RecvCounts
    eckit::mpi::Buffer<int> targetRecvCountsBuffer(comm_.size());
    comm_.allGatherv(targetRecvCounts_.begin(), targetRecvCounts_.end(), targetRecvCountsBuffer);
    std::vector<int> targetRecvCountsGlb_ = std::move(targetRecvCountsBuffer.buffer);

    // SendCounts
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      sourceSendCounts_.push_back(targetRecvCountsGlb_[jt*comm_.size()+comm_.rank()]);
    }

    // Buffer size
    sourceSendSize_ = 0;
    for (const auto & n : sourceSendCounts_) sourceSendSize_ += n;

    // SendDispls
    sourceSendDispls_.push_back(0);
    for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
      sourceSendDispls_.push_back(sourceSendDispls_[jt]+sourceSendCounts_[jt]);
    }

    // Ordered received points list
    std::vector<size_t> targetRecvOffset(comm_.size(), 0);
    std::vector<int> targetRecvPointsListOrdered(targetRecvSize_);
    for (size_t jr = 0; jr < targetRecvSize_; ++jr) {
      const size_t sourceI = targetRecvPointsList[jr]/sourceNy;
      const size_t sourceJ = targetRecvPointsList[jr]-sourceI*sourceNy;
      size_t jt = mpiTask[sourceI*sourceNy+sourceJ];
      size_t jro = targetRecvDispls_[jt]+targetRecvOffset[jt];
      targetRecvPointsListOrdered[jro] = targetRecvPointsList[jr];
      ++targetRecvOffset[jt];
    }
    std::vector<int> sourceSentPointsList(sourceSendSize_);
    comm_.allToAllv(targetRecvPointsListOrdered.data(), targetRecvCounts_.data(), targetRecvDispls_.data(),
                    sourceSentPointsList.data(), sourceSendCounts_.data(), sourceSendDispls_.data());

    // Sort indices
    std::vector<int> gij;
    for (size_t sourceJnode = 0; sourceJnode < sourceSize_; ++sourceJnode) {
      if (sourceGhostView(sourceJnode) == 0) {
        gij.push_back((sourceIndexIView(sourceJnode)-1)*sourceNy+sourceIndexJView(sourceJnode)-1);
      } else {
        gij.push_back(-1);
      }
    }
    std::vector<size_t> gidx(sourceSize_);
    std::iota(gidx.begin(), gidx.end(), 0);
    std::stable_sort(gidx.begin(), gidx.end(), [&gij](size_t i1, size_t i2)
      {return gij[i1] < gij[i2];});
    std::vector<size_t> ridx(sourceSendSize_);
    std::iota(ridx.begin(), ridx.end(), 0);
    std::stable_sort(ridx.begin(), ridx.end(), [&sourceSentPointsList](size_t i1, size_t i2)
      {return sourceSentPointsList[i1] < sourceSentPointsList[i2];});

    // Mapping for sent points
    sourceSendMapping_.resize(sourceSendSize_);
    size_t sourceJnode = 0;
    for (size_t js = 0; js < sourceSendSize_; ++js) {
      while (gij[gidx[sourceJnode]] < sourceSentPointsList[ridx[js]]) {
        ++sourceJnode;
        ASSERT(sourceJnode < sourceSize_);
      }
      sourceSendMapping_[ridx[js]] = gidx[sourceJnode];
    }

    // Sort indices
    std::vector<size_t> idx(targetRecvPointsListOrdered.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&targetRecvPointsListOrdered](size_t i1, size_t i2)
      {return targetRecvPointsListOrdered[i1] < targetRecvPointsListOrdered[i2];});

    // Compute horizontal interpolation
    for (size_t targetJnode = 0; targetJnode < targetSize_; ++targetJnode) {
      // Interpolation element default values
      std::vector<std::pair<size_t, double>> operations;

      if (targetGhostView(targetJnode) == 0) {
        // Destination grid indices
        const double targetX = targetPoints[targetJnode][0];
        bool colocatedX = false;
        int indexI = -1;
        for (size_t sourceI = 0; sourceI < sourceNx-1; ++sourceI) {
          if (std::abs(targetX-sourceGrid.x(sourceI)) < eps) {
            indexI = sourceI;
            colocatedX = true;
          }
          if (sourceGrid.x(sourceI)+eps < targetX && targetX < sourceGrid.x(sourceI+1)-eps) {
            indexI = sourceI;
            colocatedX = false;
          }
        }
        if (std::abs(targetX-sourceGrid.x(sourceNx-1)) < eps) {
          indexI = sourceNx-1;
          colocatedX = true;
        }
        const double targetY = targetPoints[targetJnode][1];
        bool colocatedY = false;
        int indexJ = -1;
        for (size_t sourceJ = 0; sourceJ < sourceNy-1; ++sourceJ) {
          if (std::abs(targetY-sourceGrid.y(sourceJ)) < eps) {
            indexJ = sourceJ;
            colocatedY = true;
          }
          if (reversedY) {
            if (sourceGrid.y(sourceJ+1)+eps < targetY && targetY < sourceGrid.y(sourceJ)-eps) {
              indexJ = sourceJ;
              colocatedY = false;
            }
          } else {
            if (sourceGrid.y(sourceJ)+eps < targetY && targetY < sourceGrid.y(sourceJ+1)-eps) {
              indexJ = sourceJ;
              colocatedY = false;
            }
          }
        }
        if (std::abs(targetY-sourceGrid.y(sourceNy-1)) < eps) {
          indexJ = sourceNy-1;
          colocatedY = true;
        }

        if (indexI == -1 || indexJ == -1) {
          // Point outside of the domain, using nearest neighbor
          if (indexI > -1) {
            if (!colocatedX &&
              (std::abs(targetX-sourceGrid.x(indexI+1)) < std::abs(targetX-sourceGrid.x(indexI)))) {
              indexI += 1;
            }
          } else {
            if (std::abs(targetX-sourceGrid.x(0)) < std::abs(targetX-sourceGrid.x(sourceNx-1))) {
              indexI = 0;
            } else {
              indexI = sourceNx-1;
            }
          }
          if (indexJ > -1) {
            if (!colocatedY &&
              (std::abs(targetY-sourceGrid.y(indexJ+1)) < std::abs(targetY-sourceGrid.y(indexJ)))) {
              indexJ += 1;
            }
          } else {
            if (std::abs(targetY-sourceGrid.y(0)) < std::abs(targetY-sourceGrid.y(sourceNy-1))) {
              indexJ = 0;
            } else {
              indexJ = sourceNy-1;
            }
            std::cout << "WARNING: point outside of the domain" << std::endl;
          }

          // Colocated point (actually nearest neighbor)
          colocatedX = true;
          colocatedY = true;
        }

        // Bilinear interpolation factor
        const double alphaX = 1.0-(sourceGrid.x(indexI)+sourceDx-targetX)/sourceDx;
        const double alphaY = reversedY ? (sourceGrid.y(indexJ)-targetY)/sourceDy
          : 1.0-(sourceGrid.y(indexJ)+sourceDy-targetY)/sourceDy;

        // Points to find
        std::vector<bool> toFind = {true, !colocatedX, !colocatedY, !colocatedX && !colocatedY};
        std::vector<size_t> valueToFind = {indexI*sourceNy+indexJ, (indexI+1)*sourceNy+indexJ,
          indexI*sourceNy+(indexJ+1), (indexI+1)*sourceNy+(indexJ+1)};
        std::vector<int> foundIndex(4, -1);

        // Binary search for each point
        for (size_t jj = 0; jj < 4; ++jj) {
          if (toFind[jj]) {
            size_t low = 0;
            size_t high = targetRecvPointsListOrdered.size()-1;
            while (low <= high) {
              size_t mid = low+(high-low)/2;
              if (valueToFind[jj] == static_cast<size_t>(targetRecvPointsListOrdered[idx[mid]])) {
                foundIndex[jj] = idx[mid];
                break;
              }
              if (valueToFind[jj] > static_cast<size_t>(targetRecvPointsListOrdered[idx[mid]])) {
                low = mid+1;
              }
              if (valueToFind[jj] < static_cast<size_t>(targetRecvPointsListOrdered[idx[mid]])) {
                high = mid-1;
              }
            }
            ASSERT(foundIndex[jj] > -1);
            ASSERT(static_cast<size_t>(targetRecvPointsListOrdered[foundIndex[jj]]) ==
              valueToFind[jj]);
          }
        }

        // Create interpolation operations
        if (colocatedX && colocatedY) {
          // Colocated point
          operations.push_back(std::make_pair(foundIndex[0], 1.0));
        } else if (colocatedY) {
          // Linear interpolation along x
          operations.push_back(std::make_pair(foundIndex[0], 1.0-alphaX));
          operations.push_back(std::make_pair(foundIndex[1], alphaX));
        } else if (colocatedX) {
          // Linear interpolation along y
          operations.push_back(std::make_pair(foundIndex[0], 1.0-alphaY));
          operations.push_back(std::make_pair(foundIndex[2], alphaY));
        } else {
          // Bilinear interpolation
          operations.push_back(std::make_pair(foundIndex[0], (1.0-alphaX)*(1.0-alphaY)));
          operations.push_back(std::make_pair(foundIndex[1], alphaX*(1.0-alphaY)));
          operations.push_back(std::make_pair(foundIndex[2], (1.0-alphaX)*alphaY));
          operations.push_back(std::make_pair(foundIndex[3], alphaX*alphaY));
        }
      }
      horInterp_.push_back(atlas::interpolation::element::InterpElement(operations));
    }
  }
}

void RegionalLinear2D::do_execute(const FieldSet& sourceFieldSet,
                                 FieldSet& targetFieldSet,
                                 Metadata& metadata) const {
  ATLAS_TRACE("atlas::interpolation::method::RegionalLinear2D::do_execute()");
  ATLAS_ASSERT(sourceFieldSet.size() == targetFieldSet.size());

  for (auto i = 0; i < sourceFieldSet.size(); ++i) {
    do_execute(sourceFieldSet[i], targetFieldSet[i], metadata);
  }
}

void RegionalLinear2D::do_execute(const Field& sourceField, Field& targetField,
                                 Metadata&) const {
  ATLAS_TRACE("atlas::interpolation::method::RegionalLinear2D::do_execute()");

  if (targetField.size() == 0) {
      return;
  }

  // Check number of levels
  ASSERT(sourceField.levels() == targetField.levels());
  const size_t nz = sourceField.levels() > 0 ? sourceField.levels() : 1;
  const size_t ndim = sourceField.levels() > 0 ? 2 : 1;

  // Scale counts and displs for all levels
  std::vector<int> sourceSendCounts3D(comm_.size());
  std::vector<int> sourceSendDispls3D(comm_.size());
  std::vector<int> targetRecvCounts3D(comm_.size());
  std::vector<int> targetRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sourceSendCounts3D[jt] = sourceSendCounts_[jt]*nz;
    sourceSendDispls3D[jt] = sourceSendDispls_[jt]*nz;
    targetRecvCounts3D[jt] = targetRecvCounts_[jt]*nz;
    targetRecvDispls3D[jt] = targetRecvDispls_[jt]*nz;
  }

  // Halo exchange
  haloExchange(sourceField);

  // Serialize
  std::vector<double> sourceSendVec(sourceSendSize_*nz);
  if (ndim == 1) {
    const auto sourceView = array::make_view<double, 1>(sourceField);
    for (size_t js = 0; js < sourceSendSize_; ++js) {
      size_t sourceJnode = sourceSendMapping_[js];
      sourceSendVec[js] = sourceView(sourceJnode);
    }
  } else if (ndim == 2) {    
    const auto sourceView = array::make_view<double, 2>(sourceField);
    for (size_t js = 0; js < sourceSendSize_; ++js) {
      for (size_t k = 0; k < nz; ++k) {
        size_t sourceJnode = sourceSendMapping_[js];
        sourceSendVec[js*nz+k] = sourceView(sourceJnode, k);
      }
    }
  }

  // Communication
  std::vector<double> targetRecvVec(targetRecvSize_*nz);
  comm_.allToAllv(sourceSendVec.data(), sourceSendCounts3D.data(), sourceSendDispls3D.data(),
                  targetRecvVec.data(), targetRecvCounts3D.data(), targetRecvDispls3D.data());

  // Interpolation
  const auto targetGhostView = array::make_view<int, 1>(targetField.functionspace().ghost());
  if (ndim == 1) {
    auto targetView = array::make_view<double, 1>(targetField);
    targetView.assign(0.0);
    for (size_t targetJnode = 0; targetJnode < targetSize_; ++targetJnode) {
      if (targetGhostView(targetJnode) == 0) {
        for (const auto & horOperation : horInterp_[targetJnode].operations()) {
          targetView(targetJnode) += horOperation.second*targetRecvVec[horOperation.first];
        }
      }
    }
  } else if (ndim == 2) {
    auto targetView = array::make_view<double, 2>(targetField);
    targetView.assign(0.0);
    for (size_t targetJnode = 0; targetJnode < targetSize_; ++targetJnode) {
      if (targetGhostView(targetJnode) == 0) {
        for (const auto & horOperation : horInterp_[targetJnode].operations()) {
          for (size_t k = 0; k < nz; ++k) {
            targetView(targetJnode, k) += horOperation.second*targetRecvVec[horOperation.first*nz+k];
          }
        }
      }
    }
  }

  // Set target field dirty
  targetField.set_dirty();
}

void RegionalLinear2D::do_execute_adjoint(FieldSet& sourceFieldSet,
                                         const FieldSet& targetFieldSet,
                                         Metadata& metadata) const {
  ATLAS_TRACE(
      "atlas::interpolation::method::RegionalLinear2D::do_execute_adjoint()");
  ATLAS_ASSERT(sourceFieldSet.size() == targetFieldSet.size());

  for (auto i = 0; i < sourceFieldSet.size(); ++i) {
    do_execute_adjoint(sourceFieldSet[i], targetFieldSet[i], metadata);
  }
}

void RegionalLinear2D::do_execute_adjoint(Field& sourceField,
                                         const Field& targetField,
                                         Metadata& metadata) const {
  ATLAS_TRACE(
      "atlas::interpolation::method::RegionalLinear2D::do_execute_adjoint()");

  if (targetField.size() == 0) {
      return;
  }

  // Check number of levels
  ASSERT(sourceField.levels() == targetField.levels());
  const size_t nz = sourceField.levels() > 0 ? sourceField.levels() : 1;
  const size_t ndim = sourceField.levels() > 0 ? 2 : 1;

  // Scale counts and displs for all levels
  std::vector<int> sourceSendCounts3D(comm_.size());
  std::vector<int> sourceSendDispls3D(comm_.size());
  std::vector<int> targetRecvCounts3D(comm_.size());
  std::vector<int> targetRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sourceSendCounts3D[jt] = sourceSendCounts_[jt]*nz;
    sourceSendDispls3D[jt] = sourceSendDispls_[jt]*nz;
    targetRecvCounts3D[jt] = targetRecvCounts_[jt]*nz;
    targetRecvDispls3D[jt] = targetRecvDispls_[jt]*nz;
  }

  // Copy destination field
  Field targetTmpField = targetField.clone();

  // Interpolation adjoint
  const auto targetGhostView = array::make_view<int, 1>(targetField.functionspace().ghost());
  std::vector<double> targetRecvVec(targetRecvSize_*nz, 0.0);
  if (ndim == 1) {
    const auto targetView = array::make_view<double, 1>(targetTmpField);
    for (size_t targetJnode = 0; targetJnode < targetSize_; ++targetJnode) {
      if (targetGhostView(targetJnode) == 0) {
        for (const auto & horOperation : horInterp_[targetJnode].operations()) {
          targetRecvVec[horOperation.first] += horOperation.second*targetView(targetJnode);
        }
      }
    }
  } else if (ndim == 2) {
    const auto targetView = array::make_view<double, 2>(targetTmpField);
    for (size_t targetJnode = 0; targetJnode < targetSize_; ++targetJnode) {
      if (targetGhostView(targetJnode) == 0) {
        for (const auto & horOperation : horInterp_[targetJnode].operations()) {
          for (size_t k = 0; k < nz; ++k) {
            targetRecvVec[horOperation.first*nz+k] += horOperation.second*targetView(targetJnode, k);
          }
        }
      }
    }
  }

  // Communication
  std::vector<double> sourceSendVec(sourceSendSize_*nz);
  comm_.allToAllv(targetRecvVec.data(), targetRecvCounts3D.data(), targetRecvDispls3D.data(),
                  sourceSendVec.data(), sourceSendCounts3D.data(), sourceSendDispls3D.data());

  // Deserialize
  if (ndim == 1) {
    auto sourceView = array::make_view<double, 1>(sourceField);
    sourceView.assign(0.0);
    for (size_t js = 0; js < sourceSendSize_; ++js) {
      size_t sourceJnode = sourceSendMapping_[js];
      sourceView(sourceJnode) += sourceSendVec[js];
    }
  } else if (ndim == 2) {
    auto sourceView = array::make_view<double, 2>(sourceField);
    sourceView.assign(0.0);
    for (size_t js = 0; js < sourceSendSize_; ++js) {
      size_t sourceJnode = sourceSendMapping_[js];
      for (size_t k = 0; k < nz; ++k) {
        sourceView(sourceJnode, k) += sourceSendVec[js*nz+k];
      }
    }
  }

  // Adjoint halo exchange
  adjointHaloExchange(sourceField);
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
