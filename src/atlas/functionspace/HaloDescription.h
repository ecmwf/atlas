#pragma once

#include <algorithm>

#include "atlas/library/config.h"

namespace atlas::functionspace {

class HaloDescription {
public:
  HaloDescription() = default;
  HaloDescription(idx_t size) :
    halo_size_(size) {
  }
  HaloDescription(idx_t begin, idx_t end) :
    halo_begin_(begin),
    halo_end_(end),
    halo_size_(halo_end_ - halo_begin_) {
  }
  template<typename View>
  HaloDescription(const View& ghost) {
    for (idx_t i=0; i<ghost.size(); ++i) {
      if (ghost[i]) {
        halo_size_ += 1;
        if (halo_begin_ == -1) {
          halo_begin_ = i;
        }
        halo_end_ = std::max(halo_end_, i+1);
      }
    }
  }

  bool contiguous() const { return halo_end_ - halo_begin_ == halo_size_; }
  idx_t size() const { return halo_size_; }
  idx_t begin() const { return halo_begin_; }
  idx_t end() const { return halo_end_; }
private:
  idx_t halo_begin_{-1};
  idx_t halo_end_{-1};
  idx_t halo_size_{0};
};

} // namespace atlas::functionspace

