#include <cuda_runtime.h>

extern "C" {
void allocate_unified_impl(double **a, int N) {
  cudaMallocManaged( a, N*sizeof(double));
}
}