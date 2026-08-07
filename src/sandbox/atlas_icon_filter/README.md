# atlas-icon-filter

This directory contains a sandbox executable for filtering an ICON NetCDF field with Atlas interpolation and spectral transforms.

## Directory Contents

- `atlas-icon-filter.F90`: Fortran executable that reads an ICON NetCDF file, interpolates `temp` to a Gaussian grid, applies a spectral cutoff filter, interpolates the filtered field back to the ICON grid, and optionally writes diagnostic outputs.
- `CMakeLists.txt`: CMake registration for the `atlas-icon-filter` executable.
- `plot_power_spectrum.py`: Helper script for plotting ASCII power-spectrum files written with `--output-spectrum`.

## Build Requirements

This executable is part of the Atlas sandbox and is only built when Atlas is configured with sandbox support:

```sh
cmake -DENABLE_SANDBOX=ON <other-options> <atlas-source-dir>
```

It is a Fortran program and therefore requires Atlas Fortran support and the mandatory Fortran-side dependencies used by this workflow:

- `fckit`
- `ectrans`
- NetCDF Fortran (`NetCDF::NetCDF_Fortran`)

The sandbox CMake file skips this executable when NetCDF Fortran is unavailable.

After configuration, build the executable with:

```sh
cmake --build <build-dir> --target atlas-icon-filter
```

## Usage

```sh
atlas-icon-filter <netcdf-file> \
                  [--output-netcdf <filename>] \
                  [--spectral-cutoff <cutoff>] \
                  [--output-spectrum] [--output-gmsh]
```

The input NetCDF file must contain ICON cell-center coordinates `clon` and `clat`, and a `temp` variable with dimensions `(ncells, plev, time)`.

Options:

- `--spectral-cutoff <cutoff>`: Total-wavenumber cutoff used by `filter_spectral_cutoff`. If omitted, the executable uses `spectral_truncation/10`.
- `--output-netcdf <filename>`: Copy the input NetCDF file and incrementally overwrite each `temp` time slice with the filtered field values. The output filename must be different from the input filename.
- `--output-spectrum`: Write ASCII power-spectrum files before and after filtering for each time step.
- `--output-gmsh`: Write Gmsh diagnostics for the ICON mesh and unfiltered/filtered fields.

## Filtering Procedure

For each time step, the executable:

1. Reads `temp` from the ICON NetCDF file into an Atlas `NodeColumns` field.
2. Interpolates from the ICON unstructured grid to a regular Gaussian grid.
3. Transforms from Gaussian grid-point space to spectral space.
4. Applies the spectral cutoff filter.
5. Transforms from spectral space back to the Gaussian grid.
6. Interpolates the filtered field back to the ICON grid.
7. Optionally writes the filtered ICON field into the copied NetCDF output.

The executable reports wall-clock timings for NetCDF reading, each filtering stage, and the combined filtering procedure for every time step.
