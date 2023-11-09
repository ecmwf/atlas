#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "eckit/mpi/Comm.h"
#include "tests/AtlasTestEnvironment.h"

#include <stdio.h>
#include <iostream>

namespace atlas {
namespace test {

CASE("test_setcomm") {
    printf(" eckit::mpi::comm () = 0x%p, eckit::mpi::comm ().size () = %8zu\n", &eckit::mpi::comm(),
           eckit::mpi::comm().size());
    std::cout << eckit::mpi::comm().name() << std::endl;
    printf(" atlas::mpi::comm () = 0x%p, atlas::mpi::comm ().size () = %8zu\n", &atlas::mpi::comm(),
           atlas::mpi::comm().size());
    std::cout << atlas::mpi::comm().name() << std::endl;

    EXPECT_EQ(eckit::mpi::comm().size(), atlas::mpi::comm().size());

    int irank = eckit::mpi::comm().rank();

    eckit::mpi::comm().split(irank, "SPLIT");

    eckit::mpi::setCommDefault("SPLIT");

    printf(" eckit::mpi::comm () = 0x%p, eckit::mpi::comm ().size () = %8zu\n", &eckit::mpi::comm(),
           eckit::mpi::comm().size());
    std::cout << eckit::mpi::comm().name() << std::endl;
    printf(" atlas::mpi::comm () = 0x%p, atlas::mpi::comm ().size () = %8zu\n", &atlas::mpi::comm(),
           atlas::mpi::comm().size());
    std::cout << atlas::mpi::comm().name() << std::endl;

    EXPECT_EQ(eckit::mpi::comm().size(), atlas::mpi::comm().size());

    std::cout << "----- STOP -----" << std::endl;
}

CASE("test_split") {
    int irank = mpi::comm("world").rank();
    auto& split_comm = mpi::comm("world").split(irank%2, "test_split");

    EXPECT_EQ(mpi::comm("world").name(), "world");
    EXPECT_EQ(split_comm.name(), "test_split");

    if( mpi::comm("world").size()%2 == 0) { // even number
       EXPECT_EQ( split_comm.size() , mpi::comm("world").size()/2);
    }

    auto& split_comm_by_name = mpi::comm("test_split");
    EXPECT_EQ( &split_comm, &split_comm_by_name);
}



}  // namespace test
}  // namespace atlas

int main(int argc, char* argv[]) {
    return atlas::test::run(argc, argv);
}
