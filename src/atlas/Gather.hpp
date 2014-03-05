#ifndef Gather_hpp
#define Gather_hpp

#include <vector>
#include <stdexcept>
#include "MPL.hpp"

namespace atlas {

class Gather
{
public:
  Gather();
  virtual ~Gather() {}

public: // methods

  void setup(const int proc[],
             const int glb_idx[],
             const int master_glb_idx[],
             const std::vector<int>& bounds,
             int par_bound );

  template <typename DATA_TYPE>
  void execute( DATA_TYPE locfield[], DATA_TYPE glbfield[], int nb_vars ) const;

  int glb_dof() const { return recvcnt_; }

private: // data
  int               packet_size_;
  int               sendcnt_;
  int               recvcnt_;
  std::vector<int>  sendcounts_;
  std::vector<int>  senddispls_;
  std::vector<int>  recvcounts_;
  std::vector<int>  recvdispls_;
  std::vector<int>  sendmap_;
  std::vector<int>  recvmap_;

  int nproc;
  int myproc;
  int root;

  std::vector<int> bounds_;
  int par_bound_;
  bool is_setup_;
};

template<typename DATA_TYPE>
void Gather::execute( DATA_TYPE locfield[], DATA_TYPE glbfield[], int nb_vars ) const
{
  if( ! is_setup_ )
  {
    throw std::runtime_error("Gather was not setup");
  }
#define FIELD_CONTIGUOUS true

  int tag=1;
  int ierr;
  int ibuf;
  int point_size = packet_size_ * nb_vars;
  int send_size = sendcnt_ * point_size;
  int recv_size = recvcnt_ * point_size;

#ifndef STACK_ARRAYS
  std::vector<DATA_TYPE> send_buffer(send_size);
  std::vector<DATA_TYPE> recv_buffer(recv_size);
  std::vector<MPI_Request> send_req(nproc);
  std::vector<MPI_Request> recv_req(nproc);
  std::vector<int> recv_displs(nproc);
  std::vector<int> recv_counts(nproc);
#else // This seems to be slower on Intel ICC 13.0.1
  DATA_TYPE send_buffer[send_size];
  DATA_TYPE recv_buffer[recv_size];
  MPI_Request send_req[nproc];
  MPI_Request recv_req[nproc];
  int recv_displs[nproc];
  int recv_counts[nproc];
#endif

  for (int jproc=0; jproc<nproc; ++jproc)
  {
    recv_counts[jproc] = recvcounts_[jproc]*point_size;
    recv_displs[jproc] = recvdispls_[jproc]*point_size;
  }

  /// -----------------------------------------------------------
  /// With mappings and everything in place, we can now call MPI

  /// Pack
  // Use original mapping + contiguous bits
  ibuf = 0;
  for( int jj=0; jj<sendcnt_; ++jj)
  {
    const int ii = point_size*sendmap_[jj];
    for( int ip=0; ip<point_size; ++ip )
      send_buffer[ibuf++] = locfield[ ii + ip ];
  }

  /// Gather
  ierr = MPI_Gatherv( send_buffer.data(), send_size, MPL::TYPE<DATA_TYPE>(),
                      recv_buffer.data(), const_cast<int*>(recv_counts.data()), const_cast<int*>(recv_displs.data()),  MPL::TYPE<DATA_TYPE>(),
                      root, MPI_COMM_WORLD );

  /// Unpack
  // Use original mapping + contiguous bits
  ibuf = 0;
  for( int jj=0; jj<recvcnt_; ++jj)
  {
    const int ii = point_size*recvmap_[jj];
    for( int ip=0; ip<point_size; ++ip)
      glbfield[ ii + ip ] = recv_buffer[ibuf++];
  }
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  Gather* atlas__Gather__new ();
  void atlas__Gather__delete (Gather* This);
  void atlas__Gather__setup (Gather* This, int proc[], int glb_idx[], int master_glb_idx[], int bounds[], int nb_bounds, int par_bound);
  void atlas__Gather__execute_int (Gather* This, int locfield[], int glbfield[], int nb_vars);
  void atlas__Gather__execute_float (Gather* This, float locfield[], float glbfield[], int nb_vars);
  void atlas__Gather__execute_double (Gather* This, double locfield[], double glbfield[], int nb_vars);
}
// ------------------------------------------------------------------


} // namespace atlas

#endif // Gather_hpp
