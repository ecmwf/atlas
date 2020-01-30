Usage with CMake   {#cmake} 
================

Enabling AddressSanitizer (for debugging)
-----------------------------------------

Tested compiler: clang-7

cmake options:

    -DBUILD_fckit=OFF -DENABLE_FORTRAN=OFF \
    -DCMAKE_CXX_FLAGS='-fsanitize=address -fno-omit-frame-pointer' \
    -DCMAKE_C_FLAGS='-fsanitize=address -fno-omit-frame-pointer' \
    -DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address

create a suppressions file (e.g. asan.supp)

    leak:libomp
    leak:libopen-pal
    leak:MPI_Init

environment variable:

    export LSAN_OPTIONS=suppressions=<path-to-asan.supp>:fast_unwind_on_malloc=0

In `LSAN_OPTIONS`, the option `fast_unwind_on_malloc=0` is required to get full stacktrace to be able to suppress MPI_Init
