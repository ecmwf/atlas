#include <string>
#include <vector>

#include "eckit/log/Bytes.h"
#include "eckit/log/ResourceUsage.h"
#include "eckit/system/ResourceUsage.h"

#include "atlas/domain.h"
#include "atlas/grid.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/trace/StopWatch.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

// Computation of gaussian latitudes for large grids can be quite costly in terms of
// time and memory.
// Atlas has a number of commonly used resolutions precomputed which can speed up Gaussian grids
// and prevents large peaks in memory usage. See gaussian Spacing.
//
// If more high resolution Gaussian latitudes need to be precomputed they can be added there, or in
// a dynamic library with self-registration.

struct Trace {
    struct Indentor {
        Indentor() { Log::info().indent(); }
        ~Indentor() { Log::info().unindent(); }
    } indentor;
    eckit::ResourceUsage resource_usage;
    atlas::Trace trace;

    Trace(const eckit::CodeLocation& where, const std::string& what):
        indentor{}, resource_usage{what, Log::debug()}, trace{where, what} {}

    double elapsed() const { return trace.elapsed(); }
    static size_t peakMemory() { return eckit::system::ResourceUsage().maxResidentSetSize(); }
    void stopAndReport() {
        trace.stop();
        Log::info() << "Time taken:        " << elapsed() << " sec" << std::endl;
        Log::info() << "Peak memory usage: " << eckit::Bytes(peakMemory()) << std::endl;
    }
};

constexpr size_t MB = 1024 * 1024;

CASE("test resources for cropping large grids") {
    std::vector<std::string> gridnames{"L40000x20000", "N8000", "O8000"};

    for (auto& gridname : gridnames) {
        EXPECT(Trace::peakMemory() < 100 * MB);

        SECTION(std::string("section") + gridname) {
            Trace trace(Here(), "Grid{" + gridname + ", GlobalDomain{-180}}");
            auto grid = Grid{gridname, GlobalDomain{-180}};
            trace.stopAndReport();
            EXPECT(trace.elapsed() < 10.);
            EXPECT(trace.peakMemory() < 100 * MB);
        }
    }
}


}  // namespace test
}  // namespace atlas


int main(int argc, char* argv[]) {
    return atlas::test::run(argc, argv);
}
