#ifndef Comm_hpp
#define Comm_hpp

#include <mpi.h>
#include <vector>

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

  int index(int i, int j, int k, int ni, int nj, int nk) const
  {
    return( i + ni*( j + nj*k) );
  }

  int index(int i, int j, int ni, int nj) const
  {
    return( i + ni*j );
  }


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
};

template<typename DATA_TYPE>
void HaloExchange::execute( DATA_TYPE field[], int nb_vars ) const
{
  using namespace detail;
  std::vector<MPI_Request> send_req(nproc);
  std::vector<MPI_Request> recv_req(nproc);
  std::vector<int> send_displs(nproc);
  std::vector<int> recv_displs(nproc);
  std::vector<int> send_counts(nproc);
  std::vector<int> recv_counts(nproc);
  int tag=1;
  int ierr;
  int send_size = sync_sendcnt_ * packet_size_ * nb_vars;
  int recv_size = sync_recvcnt_ * packet_size_ * nb_vars;
  int point_size = packet_size_ * nb_vars;

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
  std::vector<int> send_map(send_size);
  std::vector<int> recv_map(recv_size);
  create_mappings(send_map,recv_map,nb_vars);

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
  std::vector<DATA_TYPE> send_buffer(send_size);
  std::vector<DATA_TYPE> recv_buffer(recv_size);

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
  for( int ii=0; ii<send_size; ++ii)
    send_buffer[ii] = field[ send_map[ii] ];

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
  for( int ii=0; ii<recv_size; ++ii)
  {
    field[ recv_map[ii] ] = recv_buffer[ii];
  }

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
