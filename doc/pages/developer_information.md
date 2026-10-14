# Developer Information {#developer_information}

## Enabling AddressSanitizer (for Debugging)

- Tested compiler: `clang-7`

### CMake Options

```sh
-DBUILD_fckit=OFF -DENABLE_FORTRAN=OFF \
-DCMAKE_CXX_FLAGS='-fsanitize=address -fno-omit-frame-pointer' \
-DCMAKE_C_FLAGS='-fsanitize=address -fno-omit-frame-pointer' \
-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address
```

Create a suppressions file (for example as in `atlas/tools/lsan.supp`):

```text
leak:libomp
leak:libopen-pal
leak:MPI_Init
```

Environment variable:

```sh
export LSAN_OPTIONS=suppressions=<path-to-lsan.supp>:fast_unwind_on_malloc=0
```

In `LSAN_OPTIONS`, the option `fast_unwind_on_malloc=0` is required to get a full stack trace in order to suppress `MPI_Init`.

## Enabling UndefinedBehaviorSanitizer (for Debugging)

- Tested compiler: `clang-7`

### CMake Options

```sh
-DCMAKE_CXX_FLAGS='-fsanitize=undefined -fno-omit-frame-pointer -fno-sanitize-recover=all -fsanitize-blacklist=<path-to-ubsan.blacklist>'
```

The `ubsan.blacklist` file contains routines or files where no errors should be reported, for example:

```text
src:*/CGAL/Compact_container.h
```

The third-party library CGAL (version 4.9) contains a possible error that is typically suppressed this way.
