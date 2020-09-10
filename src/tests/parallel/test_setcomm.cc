#include "eckit/mpi/Comm.h"
#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/parallel/mpi/mpi.h"
#include "tests/AtlasTestEnvironment.h"

#include <iostream>
#include <stdio.h>

namespace atlas {
namespace test {

CASE ("test_setcomm")
{

  printf (" eckit::mpi::comm () = 0x%llx, eckit::mpi::comm ().size () = %8d\n", &eckit::mpi::comm (), eckit::mpi::comm ().size ());
  std::cout << eckit::mpi::comm ().name () << std::endl;
  printf (" atlas::mpi::comm () = 0x%llx, atlas::mpi::comm ().size () = %8d\n", &atlas::mpi::comm (), atlas::mpi::comm ().size ());
  std::cout << atlas::mpi::comm ().name () << std::endl;

  EXPECT_EQ (eckit::mpi::comm ().size (), atlas::mpi::comm ().size ());

  int irank = eckit::mpi::comm ().rank (); 

  auto & comm = eckit::mpi::comm ().split (irank, "SPLIT");

  eckit::mpi::setCommDefault ("SPLIT");

  printf (" eckit::mpi::comm () = 0x%llx, eckit::mpi::comm ().size () = %8d\n", &eckit::mpi::comm (), eckit::mpi::comm ().size ());
  std::cout << eckit::mpi::comm ().name () << std::endl;
  printf (" atlas::mpi::comm () = 0x%llx, atlas::mpi::comm ().size () = %8d\n", &atlas::mpi::comm (), atlas::mpi::comm ().size ());
  std::cout << atlas::mpi::comm ().name () << std::endl;

  EXPECT_EQ (eckit::mpi::comm ().size (), atlas::mpi::comm ().size ());

  std::cout << "----- STOP -----" << std::endl;

}

}
}

int main (int argc, char* argv[]) 
{
  return atlas::test::run (argc, argv);
}
