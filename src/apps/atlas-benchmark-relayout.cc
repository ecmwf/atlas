/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/**
 * @file atlas-benchmark-relayout.cc
 * @brief Benchmark Atlas host/device relayout operations between blocked and nonblocked fields.
 *
 * This executable is intended as a focused performance probe for the relayout routines declared
 * in atlas/util/relayout.h.  The relayout routines copy data between Atlas fields using two
 * horizontal memory organisations:
 *
 *   - nonblocked layout, where the horizontal point dimension is a single contiguous logical
 *     dimension, for example [npts], [npts, nlev], or [npts, nlev, nvar];
 *   - blocked layout, where the horizontal point dimension is split into blocks of length
 *     nproma, for example [nblks, nproma], [nblks, nlev, nproma], or
 *     [nblks, nvar, nlev, nproma].
 *
 * The benchmark constructs one nonblocked field, one blocked field with nproma, and one second
 * blocked field with nproma_other.  It then measures three operations:
 *
 *   - blocked_to_nonblocked: copy from the blocked field to the nonblocked field;
 *   - nonblocked_to_blocked: copy from the nonblocked field to the blocked field;
 *   - blocked_to_blocked: copy between the two blocked fields, allowing different nproma values.
 *
 * Each operation is warmed up first and then timed for the requested number of iterations.  The
 * reported timing statistics are min, max, average, and sample standard deviation.  The benchmark
 * also reports two throughput metrics:
 *
 *   - memory bandwidth in GB/s, computed from the logical amount of data read and written;
 *   - element throughput in Gelem/s, computed from the logical number of element transfers.
 *
 * The bandwidth is intentionally based on logical data movement rather than allocated blocked
 * storage.  This makes runs easier to compare when nproma changes and the final block contains
 * padding elements that are allocated but not part of the logical field.  The allocated element
 * counts are still printed so that padding overhead remains visible.  The Data section also prints
 * the memory footprint of one full blocked block for both nproma and nproma_other.  This is the
 * number of bytes in one blocked slice after the leading nblk dimension, and can be compared with
 * the reported L1/L2/L3 CPU cache sizes to judge whether a full block can fit in a given cache.
 *
 * Problem parameters
 * ------------------
 *
 *   --precision=float|double
 *       Selects the value type stored in the benchmark fields.  "single" is accepted as an alias
 *       for "float".  Comparing float and double can help distinguish byte-bandwidth limits from
 *       element-rate, indexing, vectorisation, cache, or instruction-throughput limits.
 *
 *   --npts=<integer>
 *       Number of logical horizontal points.  This is the length of the first dimension in the
 *       nonblocked field.  The number of blocks in blocked fields is ceil(npts / nproma) or
 *       ceil(npts / nproma_other).  The last block may therefore be partial.
 *
 *   --nlev=<integer>
 *       Number of vertical levels.  Use 0 to omit this dimension.  If either nlev or nvar is
 *       greater than zero while the other is zero, the benchmark creates rank-2/rank-3 fields with
 *       one extra non-horizontal dimension of length max(nlev, nvar).  If both are greater than
 *       zero, the benchmark creates rank-3/rank-4 fields with both level and variable dimensions.
 *
 *   --nvar=<integer>
 *       Number of variables.  Use 0 to omit this dimension.  Together with nlev, this controls
 *       whether rank-1/rank-2, rank-2/rank-3, or rank-3/rank-4 relayout implementations are timed.
 *
 *   --nproma=<integer>
 *       Block length for the primary blocked field.  This value controls the innermost horizontal
 *       block size used by blocked_to_nonblocked and nonblocked_to_blocked.  It is also the source
 *       blocked nproma for blocked_to_blocked.
 *
 *   --nproma-other=<integer>
 *       Block length for the second blocked field used by blocked_to_blocked.  Choosing a different
 *       value from nproma measures repacking between two blocked layouts with different block
 *       lengths.
 *
 *   --on-device=true|false
 *       Selects device execution for the public relayout API.  The benchmark synchronizes the
 *       device before and after timed regions so elapsed time includes completed device work.  The
 *       host optimisation options documented below are relevant to host paths; device paths may use
 *       different implementations.
 *
 * Benchmark parameters
 * --------------------
 *
 *   --iterations=<integer>
 *       Number of timed repetitions for each relayout operation.  The reported min, max, average,
 *       and standard deviation are computed from these samples.
 *
 *   --warmup=<integer>
 *       Number of untimed repetitions before each measured operation.  Warmup helps reduce the
 *       effect of first-use overheads such as allocation side effects, cache state, lazy runtime
 *       initialisation, and device setup.
 *
 *   --verbose=true|false, -v
 *       Prints per-iteration progress in table output mode.  Verbose progress is suppressed when
 *       --format=json is selected so that the complete program output remains valid JSON.
 *
 *   --format=table|json
 *       Selects the output format.  The default table format prints human-readable runtime,
 *       problem, benchmark, optimisation, data, and result sections.  The json format prints a
 *       single JSON document containing the same information under runtime, problem_settings,
 *       benchmark_settings, optimisation_settings, shapes, data, and results.  JSON is intended for
 *       scripts and automated benchmark collection.
 *
 * Optimisation parameters
 * -----------------------
 *
 * The options in this group are passed through to the host relayout implementation via environment
 * variables before the benchmark starts.  They are benchmark controls, but they are not exclusive
 * to this benchmark.  The same environment variables can be set by users or scripts when running
 * any Atlas application that calls the host relayout routines, even if atlas-benchmark-relayout is
 * not involved.  This is useful for testing or temporarily selecting an optimisation strategy in a
 * larger program without recompiling Atlas.
 *
 *   --loop-order=nproma_innermost|nproma_outermost
 *       Sets ATLAS_RELAYOUT_LOOP_ORDER.  This controls the order of loops in optimized host
 *       blocked/nonblocked copies for rank-3 and rank-4 field shapes.
 *
 *       Default: nproma_innermost.
 *       When this benchmark is run without --loop-order, it explicitly
 *       sets ATLAS_RELAYOUT_LOOP_ORDER=nproma_innermost before calling the relayout routines.  In
 *       applications that do not use this benchmark, leaving ATLAS_RELAYOUT_LOOP_ORDER unset has
 *       the same effect: the host relayout implementation falls back to nproma_innermost.
 *
 *       nproma_innermost keeps the block-relative horizontal index jrof as the innermost loop in
 *       the optimized kernels where possible.  This tends to favor contiguous access on the blocked
 *       side and exposes short fixed-width loops for vectorization when nproma is statically known.
 *
 *       nproma_outermost makes jrof an outer loop and moves the level loop inward.  This can be
 *       useful for testing whether the nonblocked side, cache reuse, or compiler vectorization is
 *       more favorable for a particular rank, shape, compiler, and processor.
 *
 *       Equivalent environment variable use outside this benchmark:
 *
 *           ATLAS_RELAYOUT_LOOP_ORDER=nproma_innermost
 *           ATLAS_RELAYOUT_LOOP_ORDER=nproma_outermost
 *
 *   --nproma-dispatch=static|runtime|runtime_full_blocks
 *       Sets ATLAS_RELAYOUT_NPROMA_DISPATCH.  This controls how the public host relayout wrapper
 *       selects optimized implementations for nproma.
 *
 *       Default: static.
 *       When this benchmark is run without --nproma-dispatch, it explicitly sets
 *       ATLAS_RELAYOUT_NPROMA_DISPATCH=static.  In applications that do not use this benchmark,
 *       leaving ATLAS_RELAYOUT_NPROMA_DISPATCH unset also selects static dispatch.
 *
 *       static uses explicit template dispatch for selected nproma values.  This gives the compiler
 *       compile-time constants for the block length in the most optimized paths, which can improve
 *       unrolling, vectorization, and address arithmetic.
 *
 *       runtime bypasses the explicit static nproma dispatch and uses implementations where nproma
 *       is a runtime value.  This is useful as a baseline for the cost of avoiding many specialized
 *       template instantiations.
 *
 *       runtime_full_blocks still avoids the top-level static nproma dispatch, but dispatches full
 *       blocks inside the runtime implementation for selected compile-time nproma values.  It is
 *       intended to separate the benefit of compile-time full-block loop bounds from the rest of
 *       the static dispatch machinery.
 *
 *       Equivalent environment variable use outside this benchmark:
 *
 *           ATLAS_RELAYOUT_NPROMA_DISPATCH=static
 *           ATLAS_RELAYOUT_NPROMA_DISPATCH=runtime
 *           ATLAS_RELAYOUT_NPROMA_DISPATCH=runtime_full_blocks
 *
 *   --blocked-to-blocked-use-memcpy=true|false
 *       Sets ATLAS_RELAYOUT_BLOCKED_TO_BLOCKED_USE_MEMCPY.  When enabled, host blocked-to-blocked
 *       copies use std::memcpy for whole contiguous copies and for contiguous chunks within each
 *       block repacking operation.
 *
 *       Default: true.
 *       When this benchmark is run without --blocked-to-blocked-use-memcpy, it
 *       explicitly sets ATLAS_RELAYOUT_BLOCKED_TO_BLOCKED_USE_MEMCPY=1.  In applications that do
 *       not use this benchmark, leaving ATLAS_RELAYOUT_BLOCKED_TO_BLOCKED_USE_MEMCPY unset also
 *       enables this memcpy path.  This matches the historically optimized behavior.  Turning it
 *       off forces scalar element copy loops in those code paths and is useful for measuring
 *       whether memcpy calls are helping for a given problem shape and platform.
 *
 *       Equivalent environment variable use outside this benchmark:
 *
 *           ATLAS_RELAYOUT_BLOCKED_TO_BLOCKED_USE_MEMCPY=1
 *           ATLAS_RELAYOUT_BLOCKED_TO_BLOCKED_USE_MEMCPY=0
 *
 *       The relayout implementation also accepts true/false and on/off spellings for boolean
 *       environment variables.
 *
 *   --blocked-nonblocked-use-memcpy=true|false
 *       Sets ATLAS_RELAYOUT_BLOCKED_NONBLOCKED_USE_MEMCPY.  When enabled, rank-2 host
 *       blocked/nonblocked copies may use std::memcpy for contiguous spans.
 *
 *       Default: false.
 *       When this benchmark is run without --blocked-nonblocked-use-memcpy, it
 *       explicitly sets ATLAS_RELAYOUT_BLOCKED_NONBLOCKED_USE_MEMCPY=0.  In applications that do
 *       not use this benchmark, leaving ATLAS_RELAYOUT_BLOCKED_NONBLOCKED_USE_MEMCPY unset also
 *       disables this memcpy path.  The default is false because earlier measurements did not show
 *       a clear benefit for this path.  It remains exposed to make that choice easy to re-evaluate
 *       on other compilers, CPUs, and shapes.
 *
 *       Equivalent environment variable use outside this benchmark:
 *
 *           ATLAS_RELAYOUT_BLOCKED_NONBLOCKED_USE_MEMCPY=1
 *           ATLAS_RELAYOUT_BLOCKED_NONBLOCKED_USE_MEMCPY=0
 *
 *   --use-mdspan=true|false
 *       Sets ATLAS_RELAYOUT_USE_MDSPAN.  When enabled, the public relayout wrappers convert Atlas
 *       views to mdspan before dispatching to the host or device relayout implementation.
 *
 *       Default: false.
 *       When this benchmark is run without --use-mdspan, it explicitly sets
 *       ATLAS_RELAYOUT_USE_MDSPAN=0.  In applications that do not use this benchmark, leaving
 *       ATLAS_RELAYOUT_USE_MDSPAN unset also disables mdspan dispatch.
 *
 *       This option is useful for comparing the Atlas view path against the mdspan-based path
 *       under the same problem shape, compiler, and optimisation settings.
 *
 *       Equivalent environment variable use outside this benchmark:
 *
 *           ATLAS_RELAYOUT_USE_MDSPAN=1
 *           ATLAS_RELAYOUT_USE_MDSPAN=0
 *
 *   --index-operator=true|false
 *       Sets ATLAS_RELAYOUT_INDEX_OPERATOR.  When enabled, host blocked/nonblocked relayout
 *       wrappers bypass the contiguous optimized path and unconditionally use the generic
 *       index-operator fallback implementation in both blocked-to-nonblocked and
 *       nonblocked-to-blocked directions.
 *
 *       Default: false
 *       When this benchmark is run without --index-operator, it explicitly sets
 *       ATLAS_RELAYOUT_INDEX_OPERATOR=0.  In applications that do not use this benchmark,
 *       leaving ATLAS_RELAYOUT_INDEX_OPERATOR unset also keeps the optimized path enabled.
 *
 *       Equivalent environment variable use outside this benchmark:
 *
 *           ATLAS_RELAYOUT_INDEX_OPERATOR=0
 *           ATLAS_RELAYOUT_INDEX_OPERATOR=1
 *
 * Interpreting results
 * --------------------
 *
 * Benchmark results can vary significantly with compiler, optimization flags, OpenMP settings,
 * memory bandwidth, cache size, nproma, field rank, and whether the final block is partial.  For
 * stable comparisons, keep the machine load controlled, use enough iterations to reduce noise, and
 * compare table or JSON output produced with the same problem settings.  The Runtime section records
 * the CPU model, logical CPU thread count, total system memory where available, L1/L2/L3 cache sizes
 * where they can be discovered, GPU model name, GPU total global memory, GPU multiprocessor count,
 * GPU warp size, GPU clock rate, and OpenMP thread count to make benchmark logs easier to interpret
 * later.
 */

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <ctime>
#include <vector>

#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

#include "hic/hic.h"

#include "pluto/pluto.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/library.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/relayout.h"

using namespace atlas;

namespace {

struct Settings {
    idx_t npts{1000000};
    idx_t nlev{137};
    idx_t nvar{0};
    idx_t nproma{32};
    idx_t nproma_other{64};
    bool on_device{false};
    bool verbose{false};
    std::string precision{"double"};
    std::string loop_order{"nproma_innermost"};
    std::string nproma_dispatch{"static"};
    bool blocked_to_blocked_use_memcpy{true};
    bool blocked_nonblocked_use_memcpy{false};
    bool use_mdspan{false};
    bool index_operator{false};
    std::string format{"table"};
    idx_t iterations{20};
    idx_t warmup{2};
    void validate() const {
        assert_positive("npts", npts);
        assert_positive("nproma", nproma);
        assert_positive("nproma-other", nproma_other);
        assert_non_negative("nlev", nlev);
        assert_non_negative("nvar", nvar);
        assert_non_negative("warmup", warmup);
        assert_positive("iterations", iterations);
        assert_one_of("precision", {"float", "single", "double"}, precision);
        assert_one_of("format", {"table", "json"}, format);
        assert_one_of("loop-order", {"nproma_innermost", "nproma_outermost"}, loop_order);
        assert_one_of("nproma-dispatch", {"static", "runtime", "runtime_full_blocks"}, nproma_dispatch);
    }

private:

    static void assert_positive(const std::string& name, idx_t value) {
        if (value <= 0) {
            throw_Exception(name + " must be greater than zero");
        }
    }
    static void assert_non_negative(const std::string& name, idx_t value) {
        if (value < 0) {
            throw_Exception(name + " must be greater than or equal to zero");
        }
    }
    static void assert_one_of(const std::string& name, const std::vector<std::string>& options, const std::string& value) {
        if (std::find(options.begin(), options.end(), value) == options.end()) {
            std::ostringstream out;
            out << name << " must be one of {";
            for (size_t i = 0; i < options.size(); ++i) {
                if (i > 0) {
                    out << ", ";
                }
                out << options[i];
            }
            out << "}";
            throw_Exception(out.str());
        }
    }
};

struct BenchmarkFields {
    Field blocked;
    Field blocked_other;
    Field nonblocked;
    idx_t logical_elements{0};
    idx_t blocked_elements{0};
    idx_t blocked_other_elements{0};
};

struct Measurement {
    std::string name;
    std::vector<double> seconds;
    double min_seconds{0.};
    double max_seconds{0.};
    double avg_seconds{0.};
    double stddev_seconds{0.};
    double bandwidth_gbs{0.};
    double element_rate_gs{0.};
};

struct CpuCacheInfo {
    std::string l1;
    std::string l2;
    std::string l3;
};

struct GpuInfo {
    std::string model_name;
    std::string total_global_memory;
    std::uint64_t total_global_memory_bytes{0};
    int multiprocessors{0};
    int warp_size{0};
    int clock_rate_khz{0};
};

struct CpuInfo {
    std::string cpu_processor;
    std::string system_memory;
    unsigned int cpu_logical_threads{0};
    CpuCacheInfo cpu_cache;
};

struct RuntimeInfo {
    std::string date_time;
    CpuInfo cpu;
    GpuInfo gpu;
    int openmp_threads{0};
};

struct BlockMemoryInfo {
    std::size_t blocked_bytes{0};
    std::size_t blocked_other_bytes{0};
};

void synchronize_if_needed(bool on_device) {
    if (on_device && pluto::devices()) {
        HIC_CALL(hicDeviceSynchronize());
    }
}


std::string format_gib(const std::uint64_t bytes) {
    std::ostringstream out;
    out << bytes << " bytes (" << std::fixed << std::setprecision(2) << static_cast<double>(bytes) / (1024. * 1024. * 1024.) << " GiB)";
    return out.str();
}

std::string format_kib(const std::uint64_t bytes) {
    std::ostringstream out;
    out << bytes << " bytes (" << std::fixed << std::setprecision(2) << static_cast<double>(bytes) / 1024. << " KiB)";
    return out.str();
}

std::string date_time() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time;
    localtime_r(&time, &local_time);
    std::ostringstream out;
    out << std::put_time(&local_time, "%Y-%m-%d %H:%M:%S %Z");
    return out.str();
}

GpuInfo gpu_info() {
    GpuInfo info;
    if (pluto::devices() == 0) {
        info.model_name = "none";
        info.total_global_memory = "unknown";
        return info;
    }
    [[maybe_unused]] int device = 0;
    HIC_CALL(hicGetDevice(&device));
    hicDeviceProp_t properties;
    HIC_CALL(hicGetDeviceProperties(&properties, device));
    info.model_name = properties.name;
    info.total_global_memory_bytes = static_cast<std::uint64_t>(properties.totalGlobalMem);
    info.total_global_memory = format_gib(info.total_global_memory_bytes);
    info.multiprocessors = properties.multiProcessorCount;
    info.warp_size = properties.warpSize;
    info.clock_rate_khz = properties.clockRate;
    return info;
}

std::string system_memory() {
#if defined(__APPLE__)
    std::uint64_t bytes = 0;
    std::size_t size = sizeof(bytes);
    if (sysctlbyname("hw.memsize", &bytes, &size, nullptr, 0) == 0 && bytes > 0) {
        return format_gib(bytes);
    }
#else
    std::ifstream meminfo("/proc/meminfo");
    std::string key;
    std::uint64_t value = 0;
    std::string unit;
    while (meminfo >> key >> value >> unit) {
        if (key == "MemTotal:") {
            return format_gib(value * 1024);
        }
    }
#endif
    return "unknown";
}

CpuCacheInfo cpu_cache_info() {
#if defined(__APPLE__)
    auto sysctl_size_string = [](const char* name) -> std::string {
        std::uint64_t bytes = 0;
        std::size_t size = sizeof(bytes);
        if (sysctlbyname(name, &bytes, &size, nullptr, 0) == 0 && bytes > 0) {
            return format_kib(bytes);
        }
        return "unknown";
    };
    return {sysctl_size_string("hw.l1dcachesize"), sysctl_size_string("hw.l2cachesize"), sysctl_size_string("hw.l3cachesize")};
#else
    auto  linux_cache_size = [](const std::string& level, const std::string& type) -> std::string {
        for (int index = 0; index < 8; ++index) {
            const std::string base = "/sys/devices/system/cpu/cpu0/cache/index" + std::to_string(index) + "/";
            std::ifstream level_file(base + "level");
            std::ifstream type_file(base + "type");
            std::ifstream size_file(base + "size");
            std::string cache_level;
            std::string cache_type;
            std::string cache_size;
            if (level_file >> cache_level && type_file >> cache_type && size_file >> cache_size && cache_level == level && cache_type == type) {
                return cache_size;
            }
        }
        return "unknown";
    };
    return {linux_cache_size("1", "Data"), linux_cache_size("2", "Unified"), linux_cache_size("3", "Unified")};
#endif
}

unsigned int cpu_logical_threads() {
    const unsigned int threads = std::thread::hardware_concurrency();
    return threads == 0 ? 0 : threads;
}

std::string cpu_model() {
#if defined(__APPLE__)
    std::size_t size = 0;
    if (sysctlbyname("machdep.cpu.brand_string", nullptr, &size, nullptr, 0) == 0 && size > 0) {
        std::string model(size, '\0');
        if (sysctlbyname("machdep.cpu.brand_string", model.data(), &size, nullptr, 0) == 0) {
            if (!model.empty() && model.back() == '\0') {
                model.pop_back();
            }
            if (!model.empty()) {
                return model;
            }
        }
    }
#endif

    auto trim_left = [](const std::string& value) {
        const auto first = value.find_first_not_of(" \t");
        return first == std::string::npos ? value : value.substr(first);
    };

    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line)) {
        const std::string key = "model name";
        if (line.compare(0, key.size(), key) == 0) {
            const auto separator = line.find(':');
            if (separator != std::string::npos) {
                const auto value = line.substr(separator + 1);
                return trim_left(value);
            }
        }
    }
    const char* model = std::getenv("EC_PLATFORM");
    if (model) {
        return model;
    }
    return "unknown";
}

template <typename Value>
Field make_field(const std::string& name, const std::vector<idx_t>& shape) {
    return Field(name, array::make_datatype<Value>(), array::ArrayShape(shape));
}

template <typename Value>
BenchmarkFields make_fields(const Settings& settings) {
    auto ceil_div = [] (idx_t a, idx_t b) {
        return (a + b - 1) / b;
    };

    const idx_t nblks = ceil_div(settings.npts, settings.nproma);
    const idx_t nblks_other = ceil_div(settings.npts, settings.nproma_other);

    BenchmarkFields fields;
    fields.logical_elements = settings.npts;
    fields.blocked_elements = nblks * settings.nproma;
    fields.blocked_other_elements = nblks_other * settings.nproma_other;
    if (settings.nlev > 0 && settings.nvar > 0) {
        fields.blocked       = make_field<Value>("blocked", {nblks, settings.nvar, settings.nlev, settings.nproma});
        fields.blocked_other = make_field<Value>("blocked_other", {nblks_other, settings.nvar, settings.nlev, settings.nproma_other});
        fields.nonblocked    = make_field<Value>("nonblocked", {settings.npts, settings.nlev, settings.nvar});
        fields.logical_elements *= settings.nlev * settings.nvar;
        fields.blocked_elements *= settings.nlev * settings.nvar;
        fields.blocked_other_elements *= settings.nlev * settings.nvar;
    }
    else if (settings.nlev > 0 || settings.nvar > 0) {
        const idx_t entries = std::max(settings.nlev, settings.nvar);
        fields.blocked       = make_field<Value>("blocked", {nblks, entries, settings.nproma});
        fields.blocked_other = make_field<Value>("blocked_other", {nblks_other, entries, settings.nproma_other});
        fields.nonblocked    = make_field<Value>("nonblocked", {settings.npts, entries});
        fields.logical_elements *= entries;
        fields.blocked_elements *= entries;
        fields.blocked_other_elements *= entries;
    }
    else {
        fields.blocked       = make_field<Value>("blocked", {nblks, settings.nproma});
        fields.blocked_other = make_field<Value>("blocked_other", {nblks_other, settings.nproma_other});
        fields.nonblocked    = make_field<Value>("nonblocked", {settings.npts});
    }

    if (fields.blocked.rank() != fields.nonblocked.rank() + 1) {
        throw_Exception("blocked field rank must be one greater than nonblocked field rank");
    }
    if (fields.blocked.rank() != fields.blocked_other.rank()) {
        throw_Exception("blocked fields must have the same rank");
    }
    return fields;
}

void prepare_device_fields(BenchmarkFields& fields, bool on_device) {
    if (on_device) {
        fields.blocked.syncDevice();
        fields.nonblocked.syncDevice();
        fields.blocked_other.allocateDevice();
    }
}

template <typename Operation>
Measurement measure(const std::string& name, Operation operation, idx_t iterations, bool on_device, bool verbose,
                    std::size_t bytes_moved, idx_t elements_moved) {
    synchronize_if_needed(on_device);

    Measurement measurement;
    measurement.name = name;
    measurement.seconds.reserve(static_cast<std::size_t>(iterations));

    for (idx_t iter = 0; iter < iterations; ++iter) {
        if (verbose) {
            Log::info() << "  [" << (iter + 1) << "/" << iterations << "] " << name << std::endl;
        }
        const auto start = std::chrono::steady_clock::now();
        operation();
        synchronize_if_needed(on_device);
        const auto stop = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(stop - start).count();
        measurement.seconds.push_back(seconds);
    }

    measurement.min_seconds = *std::min_element(measurement.seconds.begin(), measurement.seconds.end());
    measurement.max_seconds = *std::max_element(measurement.seconds.begin(), measurement.seconds.end());

    const double total_seconds = std::accumulate(measurement.seconds.begin(), measurement.seconds.end(), 0.0);
    measurement.avg_seconds = total_seconds / static_cast<double>(measurement.seconds.size());

    double variance = 0.0;
    if (measurement.seconds.size() > 1) {
        for (const double seconds : measurement.seconds) {
            const double diff = seconds - measurement.avg_seconds;
            variance += diff * diff;
        }
        variance /= static_cast<double>(measurement.seconds.size() - 1);
    }
    measurement.stddev_seconds = measurement.seconds.size() > 1 ? std::sqrt(variance) : 0.0;
    measurement.bandwidth_gbs = static_cast<double>(bytes_moved) / measurement.avg_seconds / 1.e9;
    measurement.element_rate_gs = static_cast<double>(elements_moved) / measurement.avg_seconds / 1.e9;
    return measurement;
}

template <typename Operation>
void warmup(Operation operation, idx_t iterations, bool on_device) {
    for (idx_t iter = 0; iter < iterations; ++iter) {
        operation();
    }
    synchronize_if_needed(on_device);
}

void print_shape(const Field& field) {
    Log::info() << "    shape: [";
    for (idx_t i = 0; i < field.rank(); ++i) {
        Log::info() << (i ? ", " : "") << field.shape(i);
    }
    Log::info() << "]" << std::endl;
}

std::size_t block_elements(const Field& field) {
    std::size_t elements = 1;
    for (idx_t i = 1; i < field.rank(); ++i) {
        elements *= static_cast<std::size_t>(field.shape(i));
    }
    return elements;
}

template <typename Value>
BlockMemoryInfo block_memory_info(const BenchmarkFields& fields) {
    return {block_elements(fields.blocked) * sizeof(Value), block_elements(fields.blocked_other) * sizeof(Value)};
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char ch : value) {
        switch (ch) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default: out << ch; break;
        }
    }
    return out.str();
}

std::string shape_json(const Field& field) {
    std::ostringstream out;
    out << "[";
    for (idx_t i = 0; i < field.rank(); ++i) {
        out << (i ? ", " : "") << field.shape(i);
    }
    out << "]";
    return out.str();
}

void print_results_table(const std::vector<Measurement>& results) {
    auto to_milliseconds = [](double seconds) {
        return seconds * 1.e3;
    };
    Log::info() << std::left
                << "  " << std::setw(28) << "operation"
                << std::right << std::setw(12) << "min [ms]"
                << std::setw(12) << "max [ms]"
                << std::setw(12) << "avg [ms]"
                << std::setw(14) << "stddev [ms]"
                << std::setw(16) << "GB/s"
                << std::setw(16) << "Gelem/s" << std::endl;
    Log::info() << "  " << std::string(110, '-') << std::endl;
    for (const auto& result : results) {
        Log::info() << std::left
                    << "  " << std::setw(28) << result.name
                    << std::right << std::fixed << std::setprecision(3)
                    << std::setw(12) << to_milliseconds(result.min_seconds)
                    << std::setw(12) << to_milliseconds(result.max_seconds)
                    << std::setw(12) << to_milliseconds(result.avg_seconds)
                    << std::setw(14) << to_milliseconds(result.stddev_seconds)
                    << std::setprecision(2)
                    << std::setw(16) << result.bandwidth_gbs
                    << std::setw(16) << result.element_rate_gs << std::endl;
    }
    Log::info() << std::defaultfloat << std::setprecision(6);
}

void print_results_json(const std::vector<Measurement>& results, const std::string& indent = "") {
    Log::info() << indent << "[" << std::endl;
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];
        Log::info() << indent << "  {" << std::endl;
        Log::info() << indent << "    \"operation\": \"" << json_escape(result.name) << "\"," << std::endl;
        Log::info() << indent << "    \"min_seconds\": " << result.min_seconds << "," << std::endl;
        Log::info() << indent << "    \"max_seconds\": " << result.max_seconds << "," << std::endl;
        Log::info() << indent << "    \"avg_seconds\": " << result.avg_seconds << "," << std::endl;
        Log::info() << indent << "    \"stddev_seconds\": " << result.stddev_seconds << "," << std::endl;
        Log::info() << indent << "    \"memory_bandwidth_gbs\": " << result.bandwidth_gbs << "," << std::endl;
        Log::info() << indent << "    \"element_throughput_gs\": " << result.element_rate_gs << std::endl;
        Log::info() << indent << "  }" << (i + 1 == results.size() ? "" : ",") << std::endl;
    }
    Log::info() << indent << "]" << std::endl;
}

void print_results_table_section(const std::vector<Measurement>& results) {
    Log::info() << "Results" << std::endl;
    print_results_table(results);
}

void print_benchmark_json(const RuntimeInfo& runtime, const Settings& settings, const BenchmarkFields& fields,
                          const std::size_t bytes_moved, const idx_t elements_moved, const BlockMemoryInfo& block_memory,
                          const std::vector<Measurement>& results) {
    Log::info() << "{" << std::endl;
    Log::info() << "  \"runtime\": {" << std::endl;
    Log::info() << "    \"date_time\": \"" << json_escape(runtime.date_time) << "\"," << std::endl;
    Log::info() << "    \"cpu\": {" << std::endl;
    Log::info() << "      \"processor\": \"" << json_escape(runtime.cpu.cpu_processor) << "\"," << std::endl;
    Log::info() << "      \"system_memory\": \"" << json_escape(runtime.cpu.system_memory) << "\"," << std::endl;
    Log::info() << "      \"logical_threads\": " << runtime.cpu.cpu_logical_threads << "," << std::endl;
    Log::info() << "      \"cache\": {" << std::endl;
    Log::info() << "        \"l1\": \"" << json_escape(runtime.cpu.cpu_cache.l1) << "\"," << std::endl;
    Log::info() << "        \"l2\": \"" << json_escape(runtime.cpu.cpu_cache.l2) << "\"," << std::endl;
    Log::info() << "        \"l3\": \"" << json_escape(runtime.cpu.cpu_cache.l3) << "\"" << std::endl;
    Log::info() << "      }" << std::endl;
    Log::info() << "    }," << std::endl;
    Log::info() << "    \"gpu\": {" << std::endl;
    Log::info() << "      \"model_name\": \"" << json_escape(runtime.gpu.model_name) << "\"," << std::endl;
    Log::info() << "      \"total_global_memory\": \"" << json_escape(runtime.gpu.total_global_memory) << "\"," << std::endl;
    Log::info() << "      \"total_global_memory_bytes\": " << runtime.gpu.total_global_memory_bytes << "," << std::endl;
    Log::info() << "      \"multiprocessors\": " << runtime.gpu.multiprocessors << "," << std::endl;
    Log::info() << "      \"warp_size\": " << runtime.gpu.warp_size << "," << std::endl;
    Log::info() << "      \"clock_rate_khz\": " << runtime.gpu.clock_rate_khz << std::endl;
    Log::info() << "    }," << std::endl;
    Log::info() << "    \"openmp_threads\": " << runtime.openmp_threads << std::endl;
    Log::info() << "  }," << std::endl;
    Log::info() << "  \"problem_settings\": {" << std::endl;
    Log::info() << "    \"npts\": " << settings.npts << "," << std::endl;
    Log::info() << "    \"nlev\": " << settings.nlev << "," << std::endl;
    Log::info() << "    \"nvar\": " << settings.nvar << "," << std::endl;
    Log::info() << "    \"nproma\": " << settings.nproma << "," << std::endl;
    Log::info() << "    \"nproma_other\": " << settings.nproma_other << "," << std::endl;
    Log::info() << "    \"on_device\": " << (settings.on_device ? "true" : "false") << std::endl;
    Log::info() << "  }," << std::endl;
    Log::info() << "  \"benchmark_settings\": {" << std::endl;
    Log::info() << "    \"verbose\": " << (settings.verbose ? "true" : "false") << "," << std::endl;
    Log::info() << "    \"precision\": \"" << json_escape(settings.precision) << "\"," << std::endl;
    Log::info() << "    \"format\": \"" << json_escape(settings.format) << "\"," << std::endl;
    Log::info() << "    \"iterations\": " << settings.iterations << "," << std::endl;
    Log::info() << "    \"warmup\": " << settings.warmup << std::endl;
    Log::info() << "  }," << std::endl;
    Log::info() << "  \"optimisation_settings\": {" << std::endl;
    Log::info() << "    \"loop_order\": \"" << json_escape(settings.loop_order) << "\"," << std::endl;
    Log::info() << "    \"nproma_dispatch\": \"" << json_escape(settings.nproma_dispatch) << "\"," << std::endl;
    Log::info() << "    \"blocked_to_blocked_use_memcpy\": " << (settings.blocked_to_blocked_use_memcpy ? "true" : "false") << "," << std::endl;
    Log::info() << "    \"blocked_nonblocked_use_memcpy\": " << (settings.blocked_nonblocked_use_memcpy ? "true" : "false") << "," << std::endl;
    Log::info() << "    \"use_mdspan\": " << (settings.use_mdspan ? "true" : "false") << "," << std::endl;
    Log::info() << "    \"index_operator\": " << settings.index_operator << std::endl;
    Log::info() << "  }," << std::endl;
    Log::info() << "  \"data\": {" << std::endl;
    Log::info() << "    \"value_type\": \"" << json_escape(settings.precision) << "\"," << std::endl;
    Log::info() << "    \"nonblocked\": {" << std::endl;
    Log::info() << "      \"shape\": " << shape_json(fields.nonblocked) << "," << std::endl;
    Log::info() << "      \"elements\": " << fields.logical_elements << "," << std::endl;
    Log::info() << "      \"memory_bytes\": " << bytes_moved << "," << std::endl;
    Log::info() << "      \"memory_gib\": " << static_cast<double>(bytes_moved)/ (1024. * 1024. * 1024.) << std::endl;
    Log::info() << "    }," << std::endl;
    Log::info() << "    \"blocked\": {" << std::endl;
    Log::info() << "      \"shape\": " << shape_json(fields.blocked) << "," << std::endl;
    Log::info() << "      \"elements\": " << fields.blocked_elements << "," << std::endl;
    Log::info() << "      \"memory_bytes\": " << bytes_moved << "," << std::endl;
    Log::info() << "      \"memory_gib\": " << static_cast<double>(bytes_moved)/ (1024. * 1024. * 1024.) << "," << std::endl;
    Log::info() << "      \"block_memory_bytes\": " << block_memory.blocked_bytes << "," << std::endl;
    Log::info() << "      \"block_memory_kib\": " << static_cast<double>(block_memory.blocked_bytes) / 1024. << std::endl;
    Log::info() << "    }," << std::endl;
    Log::info() << "    \"blocked_other\": {" << std::endl;
    Log::info() << "      \"shape\": " << shape_json(fields.blocked_other) << "," << std::endl;
    Log::info() << "      \"elements\": " << fields.blocked_other_elements << "," << std::endl;
    Log::info() << "      \"memory_bytes\": " << bytes_moved << "," << std::endl;
    Log::info() << "      \"memory_gib\": " << static_cast<double>(bytes_moved)/ (1024. * 1024. * 1024.) << "," << std::endl;
    Log::info() << "      \"block_memory_bytes\": " << block_memory.blocked_other_bytes << "," << std::endl;
    Log::info() << "      \"block_memory_kib\": " << static_cast<double>(block_memory.blocked_other_bytes) / 1024. << std::endl;
    Log::info() << "    }" << std::endl;
    Log::info() << "  }," << std::endl;
    Log::info() << "  \"results\": " << std::endl;
    print_results_json(results, "  ");
    Log::info() << "}" << std::endl;
}

template <typename Value>
int run_benchmark(const Settings& settings) {
    BenchmarkFields fields = make_fields<Value>(settings);
    prepare_device_fields(fields, settings.on_device);

    const std::size_t bytes_moved = static_cast<std::size_t>(fields.logical_elements) * sizeof(Value) * 2;
    const idx_t elements_moved = fields.logical_elements * 2;
    const RuntimeInfo runtime{date_time(), {cpu_model(), system_memory(), cpu_logical_threads(), cpu_cache_info()}, gpu_info(), atlas_omp_get_max_threads()};
    const BlockMemoryInfo block_memory = block_memory_info<Value>(fields);

    if (settings.format == "table") {
        Log::info() << "Runtime" << std::endl;
        Log::info() << "  date-time: " << runtime.date_time << std::endl;
        Log::info() << "  CPU" << std::endl;
        Log::info() << "    processor: " << runtime.cpu.cpu_processor << std::endl;
        Log::info() << "    logical threads: " << runtime.cpu.cpu_logical_threads << std::endl;
        Log::info() << "    system memory: " << runtime.cpu.system_memory << std::endl;
        Log::info() << "    L1 cache: " << runtime.cpu.cpu_cache.l1 << std::endl;
        Log::info() << "    L2 cache: " << runtime.cpu.cpu_cache.l2 << std::endl;
        Log::info() << "    L3 cache: " << runtime.cpu.cpu_cache.l3 << std::endl;
        Log::info() << "  GPU" << std::endl;
        Log::info() << "    model: " << runtime.gpu.model_name << std::endl;
        Log::info() << "    total global memory: " << runtime.gpu.total_global_memory << std::endl;
        Log::info() << "    multiprocessors: " << runtime.gpu.multiprocessors << std::endl;
        Log::info() << "    warp size: " << runtime.gpu.warp_size << std::endl;
        Log::info() << "    clock rate: " << runtime.gpu.clock_rate_khz << " kHz" << std::endl;
        Log::info() << "  OpenMP threads: " << runtime.openmp_threads << std::endl;
        Log::info() << std::endl;
        Log::info() << "Problem settings:" << std::endl;
        Log::info() << "  npts: " << settings.npts << std::endl;
        Log::info() << "  nlev: " << settings.nlev << std::endl;
        Log::info() << "  nvar: " << settings.nvar << std::endl;
        Log::info() << "  nproma: " << settings.nproma << std::endl;
        Log::info() << "  nproma_other: " << settings.nproma_other << std::endl;
        Log::info() << "  on_device: " << std::boolalpha << settings.on_device << std::endl;
        Log::info() << std::endl;
        Log::info() << "Benchmark settings:" << std::endl;
        Log::info() << "  verbose: " << std::boolalpha << settings.verbose << std::endl;
        Log::info() << "  precision: " << settings.precision << std::endl;
        Log::info() << "  format: " << settings.format << std::endl;
        Log::info() << "  iterations: " << settings.iterations << std::endl;
        Log::info() << "  warmup: " << settings.warmup << std::endl;
        Log::info() << std::endl;
        Log::info() << "Optimisation settings" << std::endl;
        Log::info() << "  loop_order: " << settings.loop_order << std::endl;
        Log::info() << "  nproma_dispatch: " << settings.nproma_dispatch << std::endl;
        Log::info() << "  blocked_to_blocked_use_memcpy: " << std::boolalpha << settings.blocked_to_blocked_use_memcpy << std::endl;
        Log::info() << "  blocked_nonblocked_use_memcpy: " << std::boolalpha << settings.blocked_nonblocked_use_memcpy << std::endl;
        Log::info() << "  use_mdspan: " << std::boolalpha << settings.use_mdspan << std::endl;
        Log::info() << "  index_operator: " << std::boolalpha << settings.index_operator << std::endl;
        Log::info() << std::noboolalpha;
        Log::info() << std::endl;
        Log::info() << "Data" << std::endl;
        Log::info() << "  value type: " << settings.precision << std::endl;
        Log::info() << "  logical element transfers: " << elements_moved << std::endl;
        Log::info() << "  nonblocked" << std::endl;
        print_shape(fields.nonblocked);
        Log::info() << "    elements: " << fields.logical_elements << std::endl;
        Log::info() << "    memory: " << format_gib(bytes_moved) << std::endl;
        Log::info() << "  blocked" << std::endl;
        print_shape(fields.blocked);
        Log::info() << "    elements: " << fields.blocked_elements << std::endl;
        Log::info() << "    memory: " << format_gib(bytes_moved) << std::endl;
        Log::info() << "    memory per full block: " << format_kib(block_memory.blocked_bytes) << std::endl;
        Log::info() << "  blocked_other" << std::endl;
        print_shape(fields.blocked_other);
        Log::info() << "    elements: " << fields.blocked_other_elements << std::endl;
        Log::info() << "    memory: " << format_gib(bytes_moved) << std::endl;
        Log::info() << "    memory per full block: " << format_kib(block_memory.blocked_other_bytes) << std::endl;
        Log::info() << std::endl;
    }

    warmup([&]() { copy_blocked_to_nonblocked(fields.blocked, fields.nonblocked, settings.on_device); }, settings.warmup,
           settings.on_device);
    auto b2n = measure("blocked_to_nonblocked",
                       [&]() { copy_blocked_to_nonblocked(fields.blocked, fields.nonblocked, settings.on_device); },
                       settings.iterations, settings.on_device, settings.verbose && settings.format == "table", bytes_moved, elements_moved);

    warmup([&]() { copy_nonblocked_to_blocked(fields.nonblocked, fields.blocked, settings.on_device); }, settings.warmup,
           settings.on_device);
    auto n2b = measure("nonblocked_to_blocked",
                       [&]() { copy_nonblocked_to_blocked(fields.nonblocked, fields.blocked, settings.on_device); },
                       settings.iterations, settings.on_device, settings.verbose && settings.format == "table", bytes_moved, elements_moved);

    warmup([&]() { copy_blocked_to_blocked(fields.blocked, fields.blocked_other, settings.on_device); }, settings.warmup,
           settings.on_device);
    auto b2b = measure("blocked_to_blocked",
                       [&]() { copy_blocked_to_blocked(fields.blocked, fields.blocked_other, settings.on_device); },
                       settings.iterations, settings.on_device, settings.verbose && settings.format == "table", bytes_moved, elements_moved);

    const std::vector<Measurement> results{b2n, n2b, b2b};
    if (settings.format == "json") {
        print_benchmark_json(runtime, settings, fields, bytes_moved, elements_moved, block_memory, results);
    }
    else {
        print_results_table_section(results);
    }

    return 0;
}

}  // namespace

class Program : public AtlasTool {
public:
    Program(int argc, char** argv): AtlasTool(argc, argv) {
        add_option(new Separator("Problem parameters"));
        add_option(new SimpleOption<std::string>("precision", "Value type: float or double. Default=double"));
        add_option(new SimpleOption<long>("npts", "Number of horizontal points. Default=1000000"));
        add_option(new SimpleOption<long>("nlev", "Number of levels. Use 0 to omit this dimension. Default=137"));
        add_option(new SimpleOption<long>("nvar", "Number of variables. Use 0 to omit this dimension. Default=0"));
        add_option(new SimpleOption<long>("nproma", "Blocked layout vector length. Default=32"));
        add_option(new SimpleOption<long>("nproma-other", "Blocked layout vector length for the other field. Default=8192"));
        add_option(new SimpleOption<bool>("on-device", "Run relayout on device. Default=false"));
        add_option(new Separator("Benchmark parameters"));
        add_option(new SimpleOption<long>("iterations", "Timed iterations. Default=20"));
        add_option(new SimpleOption<long>("warmup", "Warmup iterations. Default=2"));
        add_option(new SimpleOption<bool>("verbose,v", "Print per-iteration progress output. Default=false"));
        add_option(new SimpleOption<std::string>("format", "Result output format: table or json. Default=table"));
        add_option(new Separator("Optimisation parameters"));
        add_option(new SimpleOption<std::string>("loop-order", "Host relayout loop order: nproma_innermost or nproma_outermost. Default=nproma_innermost"));
        add_option(new SimpleOption<std::string>("nproma-dispatch", "Host relayout nproma dispatch: static, runtime, or runtime_full_blocks. Default=static"));
        add_option(new SimpleOption<bool>("blocked-to-blocked-use-memcpy", "Use memcpy for host blocked-to-blocked contiguous chunks. Default=true"));
        add_option(new SimpleOption<bool>("blocked-nonblocked-use-memcpy", "Use memcpy for host rank-2 blocked/nonblocked copies. Default=false"));
        add_option(new SimpleOption<bool>("use-mdspan", "Use mdspan dispatch in the public relayout wrappers. Default=false"));
        add_option(new SimpleOption<bool>("index-operator", "Force blocked/nonblocked relayout to use the generic index-operator fallback path. Default=false"));
    }

    std::string briefDescription() override { return "Benchmark relayout between blocked and nonblocked field layouts"; }
    std::string usage() override { return name() + " [OPTION]... [--help]"; }

    int execute(const Args& args) override {
        Settings settings;
        long value = 0;
        if (args.get("npts", value)) {
            settings.npts = static_cast<idx_t>(value);
        }
        if (args.get("nlev", value)) {
            settings.nlev = static_cast<idx_t>(value);
        }
        if (args.get("nvar", value)) {
            settings.nvar = static_cast<idx_t>(value);
        }
        if (args.get("nproma", value)) {
            settings.nproma = static_cast<idx_t>(value);
        }
        if (args.get("nproma-other", value)) {
            settings.nproma_other = static_cast<idx_t>(value);
        }
        if (args.get("iterations", value)) {
            settings.iterations = static_cast<idx_t>(value);
        }
        if (args.get("warmup", value)) {
            settings.warmup = static_cast<idx_t>(value);
        }
        args.get("on-device", settings.on_device);
        args.get("verbose", settings.verbose);
        args.get("precision", settings.precision);
        args.get("format", settings.format);
        args.get("loop-order", settings.loop_order);
        args.get("nproma-dispatch", settings.nproma_dispatch);
        args.get("blocked-to-blocked-use-memcpy", settings.blocked_to_blocked_use_memcpy);
        args.get("blocked-nonblocked-use-memcpy", settings.blocked_nonblocked_use_memcpy);
        args.get("use-mdspan", settings.use_mdspan);
        args.get("index-operator", settings.index_operator);

        settings.validate();

        constexpr int overwrite = 1;
        ::setenv("ATLAS_RELAYOUT_LOOP_ORDER", settings.loop_order.c_str(), overwrite);
        ::setenv("ATLAS_RELAYOUT_NPROMA_DISPATCH", settings.nproma_dispatch.c_str(), overwrite);
        ::setenv("ATLAS_RELAYOUT_BLOCKED_TO_BLOCKED_USE_MEMCPY", settings.blocked_to_blocked_use_memcpy ? "1" : "0", overwrite);
        ::setenv("ATLAS_RELAYOUT_BLOCKED_NONBLOCKED_USE_MEMCPY", settings.blocked_nonblocked_use_memcpy ? "1" : "0", overwrite);
        ::setenv("ATLAS_RELAYOUT_USE_MDSPAN", settings.use_mdspan ? "1" : "0", overwrite);
        ::setenv("ATLAS_RELAYOUT_INDEX_OPERATOR", settings.index_operator ? "1" : "0", overwrite);

        if (settings.precision == "float" || settings.precision == "single") {
            settings.precision = "float";
            return run_benchmark<float>(settings);
        }
        else if (settings.precision == "double") {
            settings.precision = "double";
            return run_benchmark<double>(settings);
        }
        else {
            throw_Exception("precision must be 'float' or 'double'");
        }
    }
};

int main(int argc, char** argv) {
    Program tool(argc, argv);
    return tool.start();
}