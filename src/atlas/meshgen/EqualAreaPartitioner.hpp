

//     Purpose.
//     --------
//           eq_regions provides the code to perform a high level 
//           partitioning of the surface of a sphere into regions of
//           equal area and small diameter.
//
//     Background.
//     -----------
//     This C++ implementation is ported from the MATLAB
//     "Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox" of the
//     same name developed by Paul Leopardi at the University of New South Wales. 
//     This version has been coded specifically for the case of partitioning the 
//     surface of a sphere or S^dim (where dim=2) as denoted in the original code.
//     Only a subset of the original eq_regions package has been coded to determin
//     the high level distribution of regions on a sphere, as the detailed 
//     distribution of grid points to each region is left to implentatios.
//
//     The following copyright notice for the eq_regions package is included from 
//     the original MatLab release.
//
//     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     + Release 1.10 2005-06-26                                                 +
//     +                                                                         +
//     + Copyright (c) 2004, 2005, University of New South Wales                 +
//     +                                                                         +
//     + Permission is hereby granted, free of charge, to any person obtaining   +
//     + a copy of this software and associated documentation files (the         +
//     + "Software"), to deal in the Software without restriction, including     +
//     + without limitation the rights to use, copy, modify, merge, publish,     +
//     + distribute, sublicense, and/or sell copies of the Software, and to      +
//     + permit persons to whom the Software is furnished to do so, subject to   +
//     + the following conditions:                                               +
//     +                                                                         +
//     + The above copyright notice and this permission notice shall be included +
//     + in all copies or substantial portions of the Software.                  +
//     +                                                                         +
//     + THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         +
//     + EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      +
//     + MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  +
//     + IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    +
//     + CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    +
//     + TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       +
//     + SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  +
//     +                                                                         +
//     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   

#ifndef EqualAreaPartitioner_hpp
#define EqualAreaPartitioner_hpp

#include <vector>

namespace atlas {
namespace meshgen {

void eq_caps(int N, std::vector<int>& n_regions, std::vector<double>& s_cap);
void eq_regions(int N, double xmin[], double xmax[], double ymin[], double ymax[]);


// Node struct that holds the longitude and latitude in millidegrees (integers)
// This structure is used in sorting algorithms, and uses less memory than
// if x and y were in double precision.
struct NodeInt
{
  int x, y;
  int n;
};

class EqualAreaPartitioner
{
public:
  EqualAreaPartitioner(int N);
  int partition(const double& x, const double& y) const;
  int band(const double& y) const;
  int sector(int band, const double& x) const;
  void area(int partition, int& band, int& sector) const;
  int nb_bands() const { return bands_.size(); }
  int nb_sectors(int band) const { return sectors_[band]; }
  void partition(int nb_nodes, NodeInt nodes[], int part[]) const;
  
private:
  int N_;
  std::vector<double> bands_;
  std::vector<int> sectors_;
};

} // namespace meshgen
} // namespace atlas

#endif // EqualAreaPartitioner_hpp
