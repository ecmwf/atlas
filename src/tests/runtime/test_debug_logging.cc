/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <sstream>
#include <cstdlib>
#include <thread>
#include <vector>

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

CASE("test_debug_synchronized_logging_basic") {
    // Test that ATLAS_DEBUG properly logs from multiple MPI ranks
    // without race conditions or garbled output.
    const auto& comm = mpi::comm();
    const int rank = comm.rank();
    const int size = comm.size();

    Log::info() << "[Rank " << rank << "] Starting basic debug test" << std::endl;

    // All ranks call ATLAS_DEBUG at roughly the same time
    // This creates contention in the logging system which could expose
    // race conditions if buffers are not properly synchronized
    ATLAS_DEBUG("All ranks logging simultaneously from basic test");

    // Verify each rank can log independently after the debug call
    Log::info() << "[Rank " << rank << "/" << size << "] Post-debug logging" << std::endl;

    ATLAS_DEBUG_VAR(backtrace());
    mpi::comm().barrier();
    Log::info() << "[Rank " << rank << "] Basic debug test completed" << std::endl;
}

CASE("test_debug_with_variable_logging") {
    // Test ATLAS_DEBUG_VAR with multiple ranks logging concurrently
    const auto& comm = mpi::comm();
    const int rank = comm.rank();
    const int size = comm.size();

    Log::info() << "[Rank " << rank << "] Starting variable debug test" << std::endl;

    // Each rank logs its rank and size - creates concurrent logging
    ATLAS_DEBUG_VAR(rank);
    ATLAS_DEBUG_VAR(size);

    Log::info() << "[Rank " << rank << "] Variable debug test completed" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_sync_basic") {
    // Test ATLAS_DEBUG_SYNC which is designed for barrier synchronization
    // The macro itself provides synchronization, so concurrent calls
    // from all ranks should not cause races.
    const auto& comm = mpi::comm();
    const int rank = comm.rank();

    Log::info() << "[Rank " << rank << "] Starting sync debug test" << std::endl;

    // All ranks call ATLAS_DEBUG_SYNC simultaneously
    // This barrier ensures synchronization and tests buffer flushing
    ATLAS_DEBUG_SYNC();

    Log::info() << "[Rank " << rank << "] Post-sync logging" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_sync_with_message") {
    // Test ATLAS_DEBUG_SYNC with message parameter from multiple ranks
    const auto& comm = mpi::comm();
    const int rank = comm.rank();

    Log::info() << "[Rank " << rank << "] Starting sync with message test" << std::endl;

    // All ranks call ATLAS_DEBUG_SYNC with a message
    std::ostringstream oss;
    oss << "rank_" << rank << "_message";
    ATLAS_DEBUG_SYNC(oss.str());

    Log::info() << "[Rank " << rank << "] Sync with message test completed" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_sync_with_communicator") {
    // Test ATLAS_DEBUG_SYNC with explicit communicator
    const auto& comm = mpi::comm();
    const int rank = comm.rank();

    Log::info() << "[Rank " << rank << "] Starting sync with comm test" << std::endl;

    // All ranks call ATLAS_DEBUG_SYNC with explicit communicator
    ATLAS_DEBUG_SYNC(comm);

    Log::info() << "[Rank " << rank << "] Sync with comm test completed" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_sync_with_communicator_and_message") {
    // Test ATLAS_DEBUG_SYNC with both communicator and message parameters
    // This is the most complex case and would be most prone to logging races
    const auto& comm = mpi::comm();
    const int rank = comm.rank();

    Log::info() << "[Rank " << rank << "] Starting sync with comm and message test" << std::endl;

    // All ranks call ATLAS_DEBUG_SYNC with both parameters
    std::ostringstream oss;
    oss << "concurrent_sync_rank_" << rank;
    ATLAS_DEBUG_SYNC(comm, oss.str());

    Log::info() << "[Rank " << rank << "] Sync with comm and message test completed" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_concurrent_multiple_calls") {
    // Test multiple concurrent ATLAS_DEBUG calls from all ranks
    // This creates higher contention on the logging system and tests buffer flushing
    const auto& comm = mpi::comm();
    const int rank = comm.rank();

    for (int iteration = 0; iteration < 3; ++iteration) {
        std::ostringstream oss;
        oss << "concurrent_iteration_" << iteration << "_rank_" << rank;

        // All ranks simultaneously call ATLAS_DEBUG
        ATLAS_DEBUG(oss.str());

        // Brief local work to simulate realistic MPI application behavior
        Log::info() << "[Rank " << rank << "] Iteration " << iteration << " work" << std::endl;
    }

    Log::info() << "[Rank " << rank << "] Concurrent test completed" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_split_communicator") {
    // Test ATLAS_DEBUG_SYNC with a split communicator to verify it works
    // with non-default MPI communicators
    const auto& world_comm = mpi::comm();
    const int rank = world_comm.rank();

    // Split communicator by rank parity (even/odd)
    int color = rank % 2;
    mpi::comm().split(color, "split_" + std::to_string(color));

    const auto& split_comm = mpi::comm("split_" + std::to_string(color));
    const int local_rank = split_comm.rank();
    const int local_size = split_comm.size();

    Log::info() << "[WorldRank " << rank << ", LocalRank " << local_rank << "/" << local_size
                << "] Starting split comm test" << std::endl;

    // Call ATLAS_DEBUG_SYNC with the split communicator
    ATLAS_DEBUG_SYNC(split_comm);

    Log::info() << "[WorldRank " << rank << ", LocalRank " << local_rank
                << "] Split comm test completed" << std::endl;

    mpi::comm().barrier();
}

CASE("test_debug_sync_var_basic") {
    // Test ATLAS_DEBUG_SYNC_VAR with just a variable (default communicator)
    const int rank = mpi::comm().rank();

    Log::info() << "[Rank " << rank << "] Starting sync var test" << std::endl;

    // All ranks call ATLAS_DEBUG_SYNC_VAR with their rank variable
    ATLAS_DEBUG_SYNC_VAR(rank);

    Log::info() << "[Rank " << rank << "] Sync var test completed" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_sync_var_with_communicator") {
    // Test ATLAS_DEBUG_SYNC_VAR with explicit communicator
    const auto& comm = mpi::comm();
    const int rank = comm.rank();
    const int size = comm.size();

    Log::info() << "[Rank " << rank << "] Starting sync var with comm test" << std::endl;

    // All ranks call ATLAS_DEBUG_SYNC_VAR with variable and communicator
    ATLAS_DEBUG_SYNC_VAR(comm, size);

    Log::info() << "[Rank " << rank << "] Sync var with comm test completed" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_sync_var_concurrent_multiple_calls") {
    // Test multiple concurrent ATLAS_DEBUG_SYNC_VAR calls to verify robustness
    const auto& comm = mpi::comm();
    const int rank = comm.rank();

    for (int iteration = 0; iteration < 3; ++iteration) {
        int value = rank * 10 + iteration;

        // All ranks simultaneously call ATLAS_DEBUG_SYNC_VAR
        ATLAS_DEBUG_SYNC_VAR(comm, value);

        Log::info() << "[Rank " << rank << "] Iteration " << iteration << " completed" << std::endl;
    }

    Log::info() << "[Rank " << rank << "] Concurrent sync var test completed" << std::endl;
    mpi::comm().barrier();
}

CASE("test_debug_sync_timeout_throws_when_configured") {
    const auto& comm = mpi::comm();
    if (comm.size() == 1) {
        Log::warning() << "Skipping this test which is designed for multiple ranks." << std::endl;
        return;
    }
    const std::string timeout_action = eckit::Resource<std::string>("$ATLAS_MPI_BARRIER_TIMEOUT_ACTION", "THROW");
    ATLAS_DEBUG_VAR(timeout_action);
    if (timeout_action != "THROW") {
        Log::warning() << "Skipping this test since ATLAS_MPI_BARRIER_TIMEOUT_ACTION=" << timeout_action << " is not set to THROW." << std::endl;
        return;
    }

    mpi::comm().split(0, "work_comm");  // before forcing the timeout
    const auto& work_comm = mpi::comm("work_comm");

    const int rank = comm.rank();
    const double timeout_seconds = eckit::Resource<double>("$ATLAS_MPI_BARRIER_TIMEOUT", 3.0);
    bool did_throw = false;
    if (rank > 0) {
         // Sleep longer than the barrier timeout to ensure rank 0 times out
        int sleep_ms = static_cast<int>((timeout_seconds + 1.) * 1000);
        std::this_thread::sleep_for( std::chrono::milliseconds(sleep_ms));
    }
    try {
        ATLAS_DEBUG_SYNC(work_comm, "forced timeout throw");
        // The internal barrier must be on a separate (work_comm) if we are expected to recover
    }
    catch (const eckit::Exception& exception) {
        did_throw = true;
        EXPECT(std::string(exception.what()).find("ATLAS_MPI_BARRIER_TIMEOUT") != std::string::npos);
        Log::info() << "Caught expected exception on rank " << rank << ": " << exception.what() << std::endl;
    }

    mpi::comm().barrier();
    EXPECT(did_throw == (rank == 0));
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
