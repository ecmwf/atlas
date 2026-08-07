! -----------------------------------------------------------------------------
! atlas-icon-filter
!
! Usage:
!   atlas-icon-filter <netcdf-file> \
!                     [--output-netcdf <filename>] \
!                     [--spectral-cutoff <cutoff>] \
!                     [--output-spectrum] [--output-gmsh] 
!
! Required input:
!   <netcdf-file>
!     ICON NetCDF input file containing the grid coordinates clon/clat and the
!     field variable temp with dimensions (ncells, plev, time). The executable
!     reads every time step of temp and writes filtered values only when an
!     output option is requested.
!
! Optional arguments:
!   --spectral-cutoff <cutoff>
!     Spectral total-wavenumber cutoff used by filter_spectral_cutoff. If this
!     option is omitted, the cutoff defaults to spectral_truncation/10.
!
!   --output-netcdf <filename>
!     Create a copy of the input NetCDF file, then incrementally overwrite the
!     temp slice for each processed time step with the filtered field values.
!     The output filename must be different from the input filename. Metadata,
!     coordinates, dimensions, and other variables are preserved from the copy.
!
!   --output-spectrum
!     Write ASCII power-spectrum files before and after filtering for each time
!     step: temp_<step>_unfiltered_spectrum.dat and
!     temp_<step>_filtered_spectrum.dat.
!
!   --output-gmsh
!     Write Gmsh diagnostics for the ICON mesh and unfiltered/filtered temp
!     fields: icon_mesh.msh, temp_<step>_unfiltered.msh, and
!     temp_<step>_filtered.msh.
!
! Filtering procedure for each time step:
!   1. Read temp from the ICON NetCDF file into an Atlas NodeColumns field.
!   2. Interpolate from the ICON unstructured grid to a regular Gaussian grid.
!   3. Transform from Gaussian grid-point space to spectral space.
!   4. Apply the spectral cutoff filter in spectral space.
!   5. Transform from spectral space back to the Gaussian grid.
!   6. Interpolate the filtered field back to the ICON grid.
!   7. Optionally write the filtered ICON field into the copied NetCDF output.
!
! Timing output:
!   Each time step reports wall-clock timings for NetCDF reading, each filtering
!   stage, and the combined filtering procedure.
! -----------------------------------------------------------------------------

module atlas_icon_filter_mod
  use atlas_module
  use fckit_mpi_module
  use netcdf
  implicit none
  contains

  subroutine read_grid_from_netcdf(grid, netcdf_path)
    character(len=4096), intent(in) :: netcdf_path
    type(atlas_Grid), intent(inout) :: grid

    real(kind=8), parameter :: pi = 3.14159265358979323846d0

    integer :: ncid
    integer :: ncells_dimid
    integer :: clon_varid
    integer :: clat_varid
    integer :: ncells

    real(kind=8), allocatable :: clon(:)
    real(kind=8), allocatable :: clat(:)

    logical :: file_exists
    inquire(file=trim(netcdf_path), exist=file_exists)
    if (.not. file_exists) then
        print *, 'NetCDF file does not exist: ', trim(netcdf_path)
        error stop 1
    end if


    call check_nf90(nf90_open(trim(netcdf_path), nf90_nowrite, ncid), 'nf90_open')
    call check_nf90(nf90_inq_dimid(ncid, 'ncells', ncells_dimid), 'nf90_inq_dimid(ncells)')
    call check_nf90(nf90_inquire_dimension(ncid, ncells_dimid, len=ncells), 'nf90_inquire_dimension(ncells)')

    allocate(clon(ncells))
    allocate(clat(ncells))

    call check_nf90(nf90_inq_varid(ncid, 'clon', clon_varid), 'nf90_inq_varid(clon)')
    call check_nf90(nf90_inq_varid(ncid, 'clat', clat_varid), 'nf90_inq_varid(clat)')
    call check_nf90(nf90_get_var(ncid, clon_varid, clon), 'nf90_get_var(clon)')
    call check_nf90(nf90_get_var(ncid, clat_varid, clat), 'nf90_get_var(clat)')
    call check_nf90(nf90_close(ncid), 'nf90_close')

    clon = clon * 180.0d0 / pi
    clat = clat * 180.0d0 / pi

    grid = atlas_UnstructuredGrid(clon, clat)

    deallocate(clon)
    deallocate(clat)
  end subroutine

  subroutine read_info_from_netcdf(nsteps, nlev, netcdf_path)
    integer, intent(out) :: nsteps
    integer, intent(out) :: nlev
    character(len=4096), intent(in) :: netcdf_path


    integer :: ncid
    integer :: plev_dimid
    integer :: time_dimid

    logical :: file_exists

    inquire(file=trim(netcdf_path), exist=file_exists)
    if (.not. file_exists) then
        print *, 'NetCDF file does not exist: ', trim(netcdf_path)
        error stop 1
    end if

    call check_nf90(nf90_open(trim(netcdf_path), nf90_nowrite, ncid), 'nf90_open')
    call check_nf90(nf90_inq_dimid(ncid, 'time', time_dimid), 'nf90_inq_dimid(time)')
    call check_nf90(nf90_inquire_dimension(ncid, time_dimid, len=nsteps), 'nf90_inquire_dimension(time)')
    call check_nf90(nf90_inq_dimid(ncid, 'plev', plev_dimid), 'nf90_inq_dimid(plev)')
    call check_nf90(nf90_inquire_dimension(ncid, plev_dimid, len=nlev), 'nf90_inquire_dimension(plev)')
    call check_nf90(nf90_close(ncid), 'nf90_close')
  end subroutine

  subroutine read_field_from_netcdf(field, tstep, netcdf_path)
    type(atlas_Field), intent(inout) :: field
    integer, intent(in) :: tstep
    character(len=4096), intent(in) :: netcdf_path

    integer :: ncid
    integer :: temp_varid
    integer :: ncells_dimid
    integer :: plev_dimid
    integer :: time_dimid
    integer :: ncells
    integer :: nlev
    integer :: nsteps
    integer :: start(3)
    integer :: count(3)

    real(kind=8), pointer :: field_data(:,:)

    logical :: file_exists

    inquire(file=trim(netcdf_path), exist=file_exists)
    if (.not. file_exists) then
        print *, 'NetCDF file does not exist: ', trim(netcdf_path)
        error stop 1
    end if

    call check_nf90(nf90_open(trim(netcdf_path), nf90_nowrite, ncid), 'nf90_open')
    call check_nf90(nf90_inq_dimid(ncid, 'time', time_dimid), 'nf90_inq_dimid(time)')
    call check_nf90(nf90_inquire_dimension(ncid, time_dimid, len=nsteps), 'nf90_inquire_dimension(time)')
    call check_nf90(nf90_inq_dimid(ncid, 'plev', plev_dimid), 'nf90_inq_dimid(plev)')
    call check_nf90(nf90_inquire_dimension(ncid, plev_dimid, len=nlev), 'nf90_inquire_dimension(plev)')
    call check_nf90(nf90_inq_dimid(ncid, 'ncells', ncells_dimid), 'nf90_inq_dimid(ncells)')
    call check_nf90(nf90_inquire_dimension(ncid, ncells_dimid, len=ncells), 'nf90_inquire_dimension(ncells)')
    call check_nf90(nf90_inq_varid(ncid, 'temp', temp_varid), 'nf90_inq_varid(temp)')

    if (nsteps < 1) then
      print *, 'Variable temp has no time steps'
      error stop 1
    end if

    if (field%rank() /= 2) then
      print *, 'Expected rank-2 Atlas field for temp, got rank = ', field%rank()
      error stop 1
    end if

    if (field%shape(1) /= nlev .or. field%shape(2) < ncells) then
      print *, 'Atlas field shape does not match NetCDF temp first time-slice'
      print *, 'field shape   = ', field%shape(1), field%shape(2)
      print *, 'expected shape= ', nlev, ncells
      error stop 1
    end if

    call field%data(field_data)
    start = [1, 1, tstep]
    count = [ncells, nlev, 1]
    call check_nf90(nf90_get_var(ncid, temp_varid, field_data(1:nlev,1:ncells), start=start, count=count), 'nf90_get_var(temp)')
    call check_nf90(nf90_close(ncid), 'nf90_close')
    call field%set_dirty()
  end subroutine

subroutine copy_file(source_path, target_path)
  character(len=*), intent(in) :: source_path
  character(len=*), intent(in) :: target_path

  integer, parameter :: buffer_size = 1048576

  integer :: input_unit
  integer :: output_unit
  integer :: io_status
  integer :: chunk_size
  integer(kind=8) :: file_size
  integer(kind=8) :: remaining
  character(len=1), allocatable :: buffer(:)
  logical :: file_exists
  type(fckit_mpi_comm) :: mpi

  mpi = fckit_mpi_comm()
  if (mpi%rank() /= 0) return

  if (trim(source_path) == trim(target_path)) then
    print *, 'Output NetCDF file must be different from input NetCDF file: ', trim(target_path)
    error stop 1
  end if

  inquire(file=trim(source_path), exist=file_exists, size=file_size)
  if (.not. file_exists) then
    print *, 'NetCDF file does not exist: ', trim(source_path)
    error stop 1
  end if

  open(newunit=input_unit, file=trim(source_path), access='stream', form='unformatted', status='old', action='read')
  open(newunit=output_unit, file=trim(target_path), access='stream', form='unformatted', status='replace', action='write')

  allocate(buffer(buffer_size))
  remaining = file_size
  do while (remaining > 0)
    chunk_size = int(min(remaining, int(buffer_size, kind=8)))
    read(input_unit, iostat=io_status) buffer(1:chunk_size)
    if (io_status /= 0) then
      print *, 'Failed to read input NetCDF file while copying: ', trim(source_path)
      error stop 1
    end if
    write(output_unit, iostat=io_status) buffer(1:chunk_size)
    if (io_status /= 0) then
      print *, 'Failed to write output NetCDF file while copying: ', trim(target_path)
      error stop 1
    end if
    remaining = remaining - chunk_size
  end do

  deallocate(buffer)
  close(input_unit)
  close(output_unit)
end subroutine copy_file

subroutine write_field_to_netcdf(field, tstep, netcdf_path)
  type(atlas_Field), intent(in) :: field
  integer, intent(in) :: tstep
  character(len=4096), intent(in) :: netcdf_path

  integer :: ncid
  integer :: temp_varid
  integer :: ncells_dimid
  integer :: plev_dimid
  integer :: time_dimid
  integer :: ncells
  integer :: nlev
  integer :: nsteps
  integer :: start(3)
  integer :: count(3)
  real(kind=8), pointer :: field_data(:,:)
  logical :: file_exists
  type(fckit_mpi_comm) :: mpi

  mpi = fckit_mpi_comm()
  if (mpi%rank() /= 0) return

  inquire(file=trim(netcdf_path), exist=file_exists)
  if (.not. file_exists) then
    print *, 'NetCDF output file does not exist: ', trim(netcdf_path)
    error stop 1
  end if

  call check_nf90(nf90_open(trim(netcdf_path), nf90_write, ncid), 'nf90_open(output)')
  call check_nf90(nf90_inq_dimid(ncid, 'time', time_dimid), 'nf90_inq_dimid(time)')
  call check_nf90(nf90_inquire_dimension(ncid, time_dimid, len=nsteps), 'nf90_inquire_dimension(time)')
  call check_nf90(nf90_inq_dimid(ncid, 'plev', plev_dimid), 'nf90_inq_dimid(plev)')
  call check_nf90(nf90_inquire_dimension(ncid, plev_dimid, len=nlev), 'nf90_inquire_dimension(plev)')
  call check_nf90(nf90_inq_dimid(ncid, 'ncells', ncells_dimid), 'nf90_inq_dimid(ncells)')
  call check_nf90(nf90_inquire_dimension(ncid, ncells_dimid, len=ncells), 'nf90_inquire_dimension(ncells)')
  call check_nf90(nf90_inq_varid(ncid, 'temp', temp_varid), 'nf90_inq_varid(temp)')

  if (tstep < 1 .or. tstep > nsteps) then
    print *, 'Requested time step is outside NetCDF time dimension: ', tstep
    error stop 1
  end if

  if (field%rank() /= 2) then
    print *, 'Expected rank-2 Atlas field for temp, got rank = ', field%rank()
    error stop 1
  end if

  if (field%shape(1) /= nlev .or. field%shape(2) < ncells) then
    print *, 'Atlas field shape does not match NetCDF temp output slice'
    print *, 'field shape   = ', field%shape(1), field%shape(2)
    print *, 'expected shape= ', nlev, ncells
    error stop 1
  end if

  call field%data(field_data)
  start = [1, 1, tstep]
  count = [ncells, nlev, 1]
  call check_nf90(nf90_put_var(ncid, temp_varid, field_data(1:nlev,1:ncells), start=start, count=count), 'nf90_put_var(temp)')
  call check_nf90(nf90_close(ncid), 'nf90_close(output)')
end subroutine write_field_to_netcdf

  subroutine check_nf90(status, context)
    integer, intent(in) :: status
    character(len=*), intent(in) :: context

    if (status /= nf90_noerr) then
      print *, trim(context)//': '//trim(nf90_strerror(status))
      error stop 1
    end if
  end subroutine check_nf90

  subroutine output_power_spectrum(spectral_fs, spectral_field, filename)
    type(atlas_functionspace_Spectral), intent(in) :: spectral_fs
    type(atlas_Field), intent(in) :: spectral_field
    character(len=*), intent(in) :: filename

    integer :: file_unit
    integer :: jlev
    integer :: n
    integer :: nlev
    integer :: truncation
    type(fckit_mpi_comm) :: mpi
    real(kind=8), allocatable :: spectrum(:,:)

    mpi = fckit_mpi_comm()
    truncation = spectral_fs%truncation()
    nlev = max(1, spectral_field%levels())
    allocate(spectrum(nlev, 0:truncation))
    call power_spectrum(spectral_fs, spectral_field, spectrum)

    if (mpi%rank() == 0) then
      open(newunit=file_unit, file=trim(filename), status='replace', action='write', form='formatted')
      write(file_unit, '(A)') '# Atlas spectral power spectrum'
      write(file_unit, '(A,A)') '# field: ', trim(spectral_field%name())
      write(file_unit, '(A,I0)') '# truncation: ', truncation
      write(file_unit, '(A,I0)') '# levels: ', nlev
      write(file_unit, '(A)') '# columns: wave_number power_spectrum_level_1 ... power_spectrum_level_N'
      do n = 0, truncation
        write(file_unit, '(I0)', advance='no') n
        do jlev = 1, nlev
          write(file_unit, '(1X,ES24.16)', advance='no') spectrum(jlev,n)
        end do
        write(file_unit, *)
      end do
      close(file_unit)
    end if

    deallocate(spectrum)
  end subroutine output_power_spectrum


  subroutine output_mesh(mesh, filename)
    type(atlas_Mesh), intent(in) :: mesh
    character(len=*), intent(in) :: filename
    type(atlas_Output) :: gmsh
    gmsh = atlas_output_Gmsh(filename, coordinates='xyz')
    call gmsh%write(mesh)
  end subroutine output_mesh
  subroutine output_field(field, filename)
    type(atlas_Field), intent(in) :: field
    character(len=*), intent(in) :: filename
    type(atlas_Output) :: gmsh
    gmsh = atlas_output_Gmsh(filename, coordinates='xyz')
    call gmsh%write(field)
  end subroutine output_field

  function wall_time() result(time)
    real(kind=8) :: time
    integer :: clock_count
    integer :: clock_rate

    call system_clock(clock_count, clock_rate)
    time = real(clock_count, kind=8) / real(clock_rate, kind=8)
  end function wall_time

  subroutine report_timing(step, label, elapsed)
    integer, intent(in) :: step
    character(len=*), intent(in) :: label
    real(kind=8), intent(in) :: elapsed

    write(*,'(A,I0,A,A,A,F12.6,A)') '    Timing step ', step, ': ', trim(label), ' = ', elapsed, ' s'
  end subroutine report_timing

end module atlas_icon_filter_mod

program atlas_icon_filter

  use atlas_module
  use atlas_icon_filter_mod

  implicit none

  real(kind=8), parameter :: pi = 3.14159265358979323846d0
  real(kind=8), parameter :: radius = 6731000.0d0
  character(len=4096) :: netcdf_path
  real(kind=8) :: icon_resolution
  type(atlas_Grid) :: icon_grid
  type(atlas_Mesh) :: icon_mesh
  type(atlas_Grid) :: gaussian_grid
  type(atlas_Interpolation) :: icon_interpolation
  integer(kind=8) :: gaussian_N
  integer(kind=4) :: spectral_truncation
  integer :: levels
  type(atlas_functionspace_NodeColumns) :: icon_fs
  type(atlas_functionspace_StructuredColumns) :: gaussian_fs
  type(atlas_functionspace_Spectral) :: spectral_fs
  
  type(atlas_Trans) :: trans
  type(atlas_Interpolation) :: interpolation_icon_to_gaussian
  type(atlas_Interpolation) :: interpolation_gaussian_to_icon
  type(atlas_Config) :: config

  type(atlas_Field) :: icon_field
  type(atlas_Field) :: gaussian_field
  type(atlas_Field) :: spectral_field

  character(len=4096) :: argument
  character(len=4096) :: netcdf_output_path
  logical :: output_spectrum
  logical :: output_gmsh
  logical :: output_netcdf
  integer :: iarg
  integer :: nsteps
  integer :: istep
  integer :: spectral_cutoff

  call atlas_initialize()

  output_spectrum = .false.
  output_gmsh = .false.
  output_netcdf = .false.
  netcdf_path = ''
  netcdf_output_path = ''
  spectral_cutoff = -1

  if (command_argument_count() < 1) then
    print *, 'Usage: atlas-icon-filter [--output-spectrum] [--output-gmsh] [--output-netcdf <filename>] [--spectral-cutoff <cutoff>] <netcdf-file>'
    error stop 1
  end if

iarg = 1
do while (iarg <= command_argument_count())
  call get_command_argument(iarg, argument)
  select case (trim(argument))
  case ('--output-spectrum')
    output_spectrum = .true.
  case ('--output-gmsh')
    output_gmsh = .true.
  case ('--output-netcdf')
    iarg = iarg + 1
    if (iarg > command_argument_count()) then
      print *, 'Usage: atlas-icon-filter [--output-spectrum] [--output-gmsh] [--output-netcdf <filename>] [--spectral-cutoff <cutoff>] <netcdf-file>'
      error stop 1
    end if
    call get_command_argument(iarg, netcdf_output_path)
    if (len_trim(netcdf_output_path) == 0) then
      print *, 'Usage: atlas-icon-filter [--output-spectrum] [--output-gmsh] [--output-netcdf <filename>] [--spectral-cutoff <cutoff>] <netcdf-file>'
      error stop 1
    end if
    output_netcdf = .true.
  case ('--spectral-cutoff')
    iarg = iarg + 1
    if (iarg > command_argument_count()) then
      print *, 'Usage: atlas-icon-filter [--output-spectrum] [--output-gmsh] [--output-netcdf <filename>] [--spectral-cutoff <cutoff>] <netcdf-file>'
      error stop 1
    end if
    call get_command_argument(iarg, argument)
    read(argument, *, err=100) spectral_cutoff
  case default
    if (len_trim(netcdf_path) == 0) then
      netcdf_path = argument
    else
      print *, 'Usage: atlas-icon-filter [--output-spectrum] [--output-gmsh] [--output-netcdf <filename>] [--spectral-cutoff <cutoff>] <netcdf-file>'
      error stop 1
    end if
  end select
  iarg = iarg + 1
end do

  if (len_trim(netcdf_path) == 0) then
    print *, 'Usage: atlas-icon-filter [--output-spectrum] [--output-gmsh] [--output-netcdf <filename>] [--spectral-cutoff <cutoff>] <netcdf-file>'
    error stop 1
  end if

  if (spectral_cutoff < -1) then
    print *, 'Spectral cutoff must be non-negative: ', spectral_cutoff
    error stop 1
  end if

  print *, 'netCDF file = ', trim(netcdf_path)

  call read_grid_from_netcdf(icon_grid, netcdf_path)
  call read_info_from_netcdf(nsteps, levels, netcdf_path)

  if (output_netcdf) then
    print *, 'netCDF output file = ', trim(netcdf_output_path)
    call copy_file(netcdf_path, netcdf_output_path)
  end if

  print *, 'icon grid size = ', icon_grid%size()
  icon_resolution = sqrt(4.0d0 * pi * radius*radius / icon_grid%size())
  print *, 'icon grid resolution [km] = ', icon_resolution / 1000.0d0

  gaussian_N = 10000000 / icon_resolution ! Approximate latitudes between pole and equator
  print *, 'gaussian N = ', gaussian_N
  spectral_truncation = 2 * gaussian_N - 1
  print *, 'spectral truncation = ', spectral_truncation
  if (spectral_cutoff == -1) then
    spectral_cutoff = spectral_truncation/10
  end if
  print *, 'spectral cutoff = ', spectral_cutoff

  gaussian_grid = atlas_RegularGaussianGrid(gaussian_N)
  print *, 'gaussian grid size = ', gaussian_grid%size()
  print *, 'gaussian equatorial resolution [km] = ', 2. * pi * radius / ( 4. * gaussian_N ) / 1000.0d0

  gaussian_fs = atlas_functionspace_StructuredColumns(grid=gaussian_grid, halo=1) ! halo=1 is needed for interpolation back
  spectral_fs = atlas_functionspace_Spectral(truncation=spectral_truncation)
  trans = atlas_Trans(gaussian_grid, spectral_truncation)

  icon_mesh = atlas_Mesh(icon_grid)

  if (output_gmsh) then
    call output_mesh(icon_mesh, 'icon_mesh.msh')
  endif

  icon_fs = atlas_functionspace_NodeColumns(mesh=icon_mesh)
  config = atlas_Config()
  call config%set('type', 'finite-element')
  interpolation_icon_to_gaussian = atlas_Interpolation(config, source=icon_fs, target=gaussian_fs)
  call config%set('type', 'structured-bilinear')
  interpolation_gaussian_to_icon = atlas_Interpolation(config, source=gaussian_fs, target=icon_fs)

  icon_field = icon_fs%create_field(name='temp', kind=atlas_real(8), levels=levels)
  gaussian_field = gaussian_fs%create_field(name='temp', kind=atlas_real(8), levels=levels)
  spectral_field = spectral_fs%create_field(name='temp', kind=atlas_real(8), levels=levels)

  do istep=1,nsteps; block
    character(len=4) :: istep_str
    real(kind=8) :: filtering_elapsed
    real(kind=8) :: filtering_start
    real(kind=8) :: timer_start
    real(kind=8) :: timer_elapsed
    write(istep_str, '(I0)') istep
    write(*,'(A,I0,A,I0)') 'Processing time step ', istep, ' of ', nsteps

    timer_start = wall_time()
    call read_field_from_netcdf(icon_field, istep, netcdf_path)
    timer_elapsed = wall_time() - timer_start
    call report_timing(istep, 'read_field_from_netcdf', timer_elapsed)

    if (output_gmsh) then
      call output_field(icon_field, 'temp_'//trim(istep_str)//'_unfiltered.msh')
    endif

    ! ---- Filtering procedure:
    ! 1. Interpolate from ICON grid to Gaussian grid
    ! 2. Transform from Gaussian grid to spectral space
    ! 3. Apply spectral cutoff filter
    ! 4. Transform back from spectral space to Gaussian grid
    ! 5. Interpolate back from Gaussian grid to ICON grid

    filtering_elapsed = 0

    timer_start = wall_time()
    call interpolation_icon_to_gaussian%execute(icon_field, gaussian_field)
    timer_elapsed = wall_time() - timer_start
    filtering_elapsed = filtering_elapsed + timer_elapsed
    call report_timing(istep, 'interpolation_icon_to_gaussian', timer_elapsed)

    timer_start = wall_time()
    call trans%dirtrans(gaussian_field, spectral_field)
    timer_elapsed = wall_time() - timer_start
    filtering_elapsed = filtering_elapsed + timer_elapsed
    call report_timing(istep, 'dirtrans', timer_elapsed)

    if (output_spectrum) then
      call output_power_spectrum(spectral_fs, spectral_field, 'temp_'//trim(istep_str)//'_unfiltered_spectrum.dat')
    end if

    timer_start = wall_time()
    call filter_spectral_cutoff(spectral_fs, spectral_field, cutoff=spectral_cutoff)
    timer_elapsed = wall_time() - timer_start
    filtering_elapsed = filtering_elapsed + timer_elapsed
    call report_timing(istep, 'filter_spectral_cutoff', timer_elapsed)

    if (output_spectrum) then
      call output_power_spectrum(spectral_fs, spectral_field, 'temp_'//trim(istep_str)//'_filtered_spectrum.dat')
    end if

    timer_start = wall_time()
    call trans%invtrans(spectral_field, gaussian_field)
    timer_elapsed = wall_time() - timer_start
    filtering_elapsed = filtering_elapsed + timer_elapsed
    call report_timing(istep, 'invtrans', timer_elapsed)

    timer_start = wall_time()
    call interpolation_gaussian_to_icon%execute(gaussian_field, icon_field)
    timer_elapsed = wall_time() - timer_start
    filtering_elapsed = filtering_elapsed + timer_elapsed
    call report_timing(istep, 'interpolation_gaussian_to_icon', timer_elapsed)

    call report_timing(istep, 'Total filtering procedure', filtering_elapsed)

    if (output_netcdf) then
      call write_field_to_netcdf(icon_field, istep, netcdf_output_path)
    end if

    if (output_gmsh) then
      call output_field(icon_field, 'temp_'//trim(istep_str)//'_filtered.msh')
    endif
  end block; end do

  call atlas_finalize()
  stop

100 print *, 'Invalid value for --spectral-cutoff: ', trim(argument)
  error stop 1

end program atlas_icon_filter
