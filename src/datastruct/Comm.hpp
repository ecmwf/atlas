#ifndef Comm_hpp
#define Comm_hpp

#include <mpi.h>
#include <vector>
#include <stdexcept>

namespace ecmwf {

namespace detail {
  template<typename DATA_TYPE>
  MPI_Datatype MPI_TYPE();
  template<> inline MPI_Datatype MPI_TYPE<int>()    { return MPI_INT; }
  template<> inline MPI_Datatype MPI_TYPE<float>()  { return MPI_FLOAT; }
  template<> inline MPI_Datatype MPI_TYPE<double>() { return MPI_DOUBLE; }
}

class HaloExchange
{
public:
  HaloExchange();
  virtual ~HaloExchange() {}

public: // methods

  void setup( const int proc[], 
              const int glb_idx[], 
              const std::vector<int>& bounds, 
              int par_bound );

  template <typename DATA_TYPE>
  void execute( DATA_TYPE field[], int nb_vars ) const;

private: // methods

  void create_mappings( std::vector<int>& send_map, std::vector<int>& recv_map, int nb_vars ) const;

  template<int N, int P>
  void create_mappings_impl( std::vector<int>& send_map, std::vector<int>& recv_map, int nb_vars ) const;

  int index(int i, int j, int k, int ni, int nj, int nk) const { return( i + ni*( j + nj*k) ); }

  int index(int i, int j, int ni, int nj) const { return( i + ni*j ); }


private: // data
  int               packet_size_;
  int               sync_sendcnt_;
  int               sync_recvcnt_;
  std::vector<int>  sync_sendcounts_;
  std::vector<int>  sync_senddispls_;
  std::vector<int>  sync_recvcounts_;
  std::vector<int>  sync_recvdispls_;
  std::vector<int>  sync_sendmap_;
  std::vector<int>  sync_recvmap_;

  int nproc;
  int myproc;

  std::vector<int> bounds_;
  int par_bound_;
  bool is_setup_;
};

template<typename DATA_TYPE>
void HaloExchange::execute( DATA_TYPE field[], int nb_vars ) const
{
  if( ! is_setup_ )
  {
    throw std::runtime_error("HaloExchange was not setup");
  }
#define FIELD_CONTIGUOUS true

  using namespace detail;
  int tag=1;
  int ierr;
  int ibuf;
  int point_size = packet_size_ * nb_vars;
  int send_size = sync_sendcnt_ * point_size;
  int recv_size = sync_recvcnt_ * point_size;

#ifndef STACK_ARRAYS
  std::vector<DATA_TYPE> send_buffer(send_size);
  std::vector<DATA_TYPE> recv_buffer(recv_size);
  std::vector<MPI_Request> send_req(nproc);
  std::vector<MPI_Request> recv_req(nproc);
  std::vector<int> send_displs(nproc);
  std::vector<int> recv_displs(nproc);
  std::vector<int> send_counts(nproc);
  std::vector<int> recv_counts(nproc);
#else // This seems to be slower on Intel ICC 13.0.1
  DATA_TYPE send_buffer[send_size];
  DATA_TYPE recv_buffer[recv_size];
  MPI_Request send_req[nproc];
  MPI_Request recv_req[nproc];
  int send_displs[nproc];
  int recv_displs[nproc];
  int send_counts[nproc];
  int recv_counts[nproc];
#endif


  // std::cout << myproc << "  :  field before = ";
  // for( int i=0; i< nb_vars*bounds_[par_bound_]; ++i)
  //   std::cout << field[i] << " ";
  // std::cout << std::endl;

  for (int jproc=0; jproc<nproc; ++jproc)
  {
    send_counts[jproc] = sync_sendcounts_[jproc]*point_size;
    recv_counts[jproc] = sync_recvcounts_[jproc]*point_size;
    send_displs[jproc] = sync_senddispls_[jproc]*point_size;
    recv_displs[jproc] = sync_recvdispls_[jproc]*point_size;
  }

#ifndef FIELD_CONTIGUOUS
  // Create additional mapping
  std::vector<int> send_map(send_size);
  std::vector<int> recv_map(recv_size);
  create_mappings(send_map,recv_map,nb_vars);
#endif
  // std::cout << myproc << "  :  send_map  = ";
  // for( int i=0; i< send_map.size(); ++i)
  //   std::cout << send_map[i] << " ";
  // std::cout << std::endl;

  // std::cout << myproc << "  :  recv_map  = ";
  // for( int i=0; i< recv_map.size(); ++i)
  //   std::cout << recv_map[i] << " ";
  // std::cout << std::endl;

  /// -----------------------------------------------------------
  /// With mappings and everything in place, we can now call MPI

  /// Let MPI know what we like to receive
  for( int jproc=0; jproc<nproc; ++jproc )
  {
    if(recv_counts[jproc] > 0)
    {
      ierr = MPI_Irecv( &recv_buffer[recv_displs[jproc]] , recv_counts[jproc],
                        MPI_TYPE<DATA_TYPE>(), jproc, tag, MPI_COMM_WORLD, &recv_req[jproc] );
    }
  }

  /// Pack
#ifndef FIELD_CONTIGUOUS
  // Use additional mapping
  for( int jj=0; jj<send_size; ++jj)
    send_buffer[jj] = field[ send_map[jj] ];
#else
  // Use original mapping + contiguous bits
  ibuf = 0;
  for( int jj=0; jj<sync_sendcnt_; ++jj)
  {
    const int ii = point_size*sync_sendmap_[jj];
    for( int ip=0; ip<point_size; ++ip )
      send_buffer[ibuf++] = field[ ii + ip ];
  }
#endif

  /// Send
  for( int jproc=0; jproc<nproc; ++jproc )
  {
    if(send_counts[jproc] > 0)
    {
      ierr = MPI_Isend( &send_buffer[send_displs[jproc]], send_counts[jproc],
                        MPI_TYPE<DATA_TYPE>(), jproc, tag, MPI_COMM_WORLD, &send_req[jproc] );
    }
  }

  /// Wait for receiving to finish
  for (int jproc=0; jproc<nproc; ++jproc)
  {
    if( sync_recvcounts_[jproc] > 0)
    {
      ierr = MPI_Wait(&recv_req[jproc], MPI_STATUS_IGNORE );
    }
  }

  /// Unpack
#ifndef FIELD_CONTIGUOUS
  // Use additional mapping
  for( int jj=0; jj<recv_size; ++jj)
  {
    field[ recv_map[jj] ] = recv_buffer[jj];
  }
#else
  // Use original mapping + contiguous bits
  ibuf = 0;
  for( int jj=0; jj<sync_recvcnt_; ++jj)
  {
    const int ii = point_size*sync_recvmap_[jj];
    for( int ip=0; ip<point_size; ++ip)
      field[ ii + ip ] = recv_buffer[ibuf++];
  }
#endif

  /// Wait for sending to finish
  for (int jproc=0; jproc<nproc; ++jproc)
  {
    if( sync_sendcounts_[jproc] > 0)
    {
      ierr = MPI_Wait(&send_req[jproc], MPI_STATUS_IGNORE );
    }
  }

  // std::cout << myproc << "  :  field after  = ";
  // for( int i=0; i< nb_vars*bounds_[par_bound_]; ++i)
  //   std::cout << field[i] << " ";
  // std::cout << std::endl;
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  HaloExchange* ecmwf__HaloExchange__new ();
  void ecmwf__HaloExchange__delete (HaloExchange* This);
  void ecmwf__HaloExchange__setup (HaloExchange* This, int proc[], int glb_idx[], int bounds[], int nb_bounds, int par_bound);
  void ecmwf__HaloExchange__execute_int (HaloExchange* This, int field[], int nb_vars);
  void ecmwf__HaloExchange__execute_float (HaloExchange* This, float field[], int nb_vars);
  void ecmwf__HaloExchange__execute_double (HaloExchange* This, double field[], int nb_vars);

}
// ------------------------------------------------------------------


} // namespace ecmwf

#endif // Comm_hpp
