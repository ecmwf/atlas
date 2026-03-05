#include "atlas/parallel/mpi/mpi.h"
#include "eckit/mpi/Comm.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas::test {

CASE("test_mpi_scope") {
    mpi::comm().split(mpi::comm().rank(), "comm1");
    mpi::comm().split(mpi::comm().rank(), "comm2");

    eckit::mpi::setCommDefault("comm1");
    EXPECT_EQ( mpi::comm().name(), "comm1" );
    {
        // Use manual push/pop style scope management
        mpi::scope::push("self");
        EXPECT_EQ( mpi::comm().name(), "self" );
        mpi::scope::pop();
    }
    EXPECT_EQ( mpi::comm().name(), "comm1" );

    eckit::mpi::setCommDefault("comm2");
    EXPECT_EQ( mpi::comm().name(), "comm2" );
    {
        // Use RAII style scope management using atlas::mpi::Scope
        mpi::Scope scope("self");
        EXPECT_EQ( mpi::comm().name(), "self" );
    }
    EXPECT_EQ( mpi::comm().name(), "comm2" );


    eckit::mpi::setCommDefault("comm1");
    EXPECT_EQ( mpi::comm().name(), "comm1" );
    {
        // Use RAII style scope management using atlas::mpi::Scope, now without explicitly specifying the comm name
        // We can edit at will manually the default comm within the scope, and it will be restored to "comm1" at the end of the scope
        mpi::Scope scope;
        EXPECT_EQ( mpi::comm().name(), "comm1" );
        eckit::mpi::setCommDefault("self");
        EXPECT_EQ( mpi::comm().name(), "self" );
    }
    EXPECT_EQ( mpi::comm().name(), "comm1" );

    {
        // Use RAII style scope management using atlas::mpi::Scope with Comm argument
        mpi::Scope scope(mpi::comm("self"));
        EXPECT_EQ( mpi::comm().name(), "self" );
    }
    EXPECT_EQ( mpi::comm().name(), "comm1" );

    eckit::mpi::setCommDefault("world");
}

}  // namespace atlas::test

#if ATLAS_HAVE_MPI
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
namespace atlas::test {
CASE ("test_mpi_scope_with_int and deep nesting") {
    constexpr bool eckit_unregister_comm_not_available = not ATLAS_ECKIT_VERSION_AT_LEAST(2, 0, 0);
    MPI_Comm comm1, comm2;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_split(MPI_COMM_WORLD, rank%2, rank, &comm1);
    MPI_Comm_split(MPI_COMM_WORLD, rank%2, rank, &comm2);
    int comm1_int = MPI_Comm_c2f(comm1);
    int comm2_int = MPI_Comm_c2f(comm2);
    int world_int = MPI_Comm_c2f(MPI_COMM_WORLD);
    EXPECT_EQ( mpi::comm().name(), "world");
    {
        mpi::Scope scope(comm1_int);
        EXPECT_EQ( mpi::comm().name(), "int." + std::to_string(comm1_int) );
        {
            mpi::Scope scope2(comm2_int);
            EXPECT_EQ( mpi::comm().name(), "int." + std::to_string(comm2_int) );
            {
                mpi::Scope scope3(comm1_int);
                EXPECT_EQ( mpi::comm().name(), "int." + std::to_string(comm1_int) );
                {
                    mpi::Scope scope4(world_int);
                    EXPECT_EQ( mpi::comm().name(), "world");
                }
                EXPECT_EQ( mpi::comm().name(), "int." + std::to_string(comm1_int) );
            }
            EXPECT( mpi::has_comm("int." + std::to_string(comm1_int)) );
        }
        EXPECT( not mpi::has_comm("int." + std::to_string(comm2_int)) || eckit_unregister_comm_not_available );
        EXPECT_EQ( mpi::comm().name(), "int." + std::to_string(comm1_int) );
    }
    EXPECT( not mpi::has_comm("int." + std::to_string(comm2_int)) || eckit_unregister_comm_not_available );
    EXPECT_EQ( mpi::comm().name(), "world");
    MPI_Comm_free(&comm1);
    MPI_Comm_free(&comm2);
}
}  // namespace atlas::test
#endif


int main(int argc, char* argv[]) {
    return atlas::test::run(argc, argv);
}
