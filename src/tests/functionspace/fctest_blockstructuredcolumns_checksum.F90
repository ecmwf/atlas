! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for BlockStructuredColumns Checksum functionality.
! Tests verify that checksums are computed correctly for fields with various
! configurations, including multi-field fieldsets and fields wrapping slices of data.
!
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_BlockStructuredColumns_checksum_fxt
use atlas_module
use atlas_functionspace_blockstructuredcolumns_module
use, intrinsic :: iso_c_binding
implicit none
character(len=1024) :: msg

interface fill_field_with_data
  module procedure fill_field_with_data_field
  module procedure fill_field_with_data_real64_2d
  module procedure fill_field_with_data_real64_3d
  module procedure fill_field_with_data_real64_4d
  module procedure fill_field_with_data_real32_2d
  module procedure fill_field_with_data_real32_3d
  module procedure fill_field_with_data_real32_4d
end interface

contains

  ! Helper to fill a field deterministically based on functionspace global_index.
  ! This makes values reproducible for different MPI partitionings.
  subroutine fill_field_with_data_field(fs, field, seed)
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: fs
    type(atlas_Field), intent(inout) :: field
    integer, intent(in) :: seed
    real(c_double), pointer :: data_2d(:,:)
    real(c_double), pointer :: data_3d(:,:,:)
    real(c_float), pointer :: data_2d_f(:,:)
    real(c_float), pointer :: data_3d_f(:,:,:)
    
    if( field%rank() == 2 ) then
      if( field%kind() == atlas_real(c_double) ) then
        call field%data(data_2d)
        call fill_field_with_data_real64_2d(fs, data_2d, seed)
      else if( field%kind() == atlas_real(c_float) ) then
        call field%data(data_2d_f)
        call fill_field_with_data_real32_2d(fs, data_2d_f, seed)
      endif
    else if( field%rank() == 3 ) then
      if( field%kind() == atlas_real(c_double) ) then
        call field%data(data_3d)
        call fill_field_with_data_real64_3d(fs, data_3d, seed)
      else if( field%kind() == atlas_real(c_float) ) then
        call field%data(data_3d_f)
        call fill_field_with_data_real32_3d(fs, data_3d_f, seed)
      endif
    endif
  end subroutine fill_field_with_data_field

  subroutine fill_field_with_data_real64_2d(fs, data, seed)
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: fs
    real(c_double), intent(inout) :: data(:,:)
    integer, intent(in) :: seed
    type(atlas_Field) :: field_global_index
    integer(ATLAS_KIND_GIDX), pointer :: global_index(:)
    integer(ATLAS_KIND_IDX) :: jblk, jpt, i

    field_global_index = fs%global_index()
    call field_global_index%data(global_index)

    data = real(seed, c_double)
    do jblk = 1, fs%nblks()
      do i = 1, fs%block_size(jblk)
        jpt = fs%block_begin(jblk) + i - 1
        data(i,jblk) = real(seed, c_double) + real(global_index(jpt), c_double)
      enddo
    enddo
    call field_global_index%final()
  end subroutine fill_field_with_data_real64_2d

  subroutine fill_field_with_data_real32_2d(fs, data, seed)
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: fs
    real(c_float), intent(inout) :: data(:,:)
    integer, intent(in) :: seed
    type(atlas_Field) :: field_global_index
    integer(ATLAS_KIND_GIDX), pointer :: global_index(:)
    integer(ATLAS_KIND_IDX) :: jblk, jpt, i

    field_global_index = fs%global_index()
    call field_global_index%data(global_index)

    data = real(seed, c_float)
    do jblk = 1, fs%nblks()
      do i = 1, fs%block_size(jblk)
        jpt = fs%block_begin(jblk) + i - 1
        data(i,jblk) = real(seed, c_float) + real(global_index(jpt), c_float)
      enddo
    enddo
    call field_global_index%final()
  end subroutine fill_field_with_data_real32_2d

  subroutine fill_field_with_data_real64_3d(fs, data, seed)
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: fs
    real(c_double), intent(inout) :: data(:,:,:)
    integer, intent(in) :: seed
    type(atlas_Field) :: field_global_index
    integer(ATLAS_KIND_GIDX), pointer :: global_index(:)
    integer(ATLAS_KIND_IDX) :: jblk, jpt, i
    integer :: j

    field_global_index = fs%global_index()
    call field_global_index%data(global_index)

    data = real(seed, c_double)
    do jblk = 1, fs%nblks()
      do i = 1, fs%block_size(jblk)
        jpt = fs%block_begin(jblk) + i - 1
        do j = 1, size(data,2)
          data(i,j,jblk) = real(seed, c_double) + real(global_index(jpt), c_double) + real(j, c_double)
        enddo
      enddo
    enddo
    call field_global_index%final()
  end subroutine fill_field_with_data_real64_3d

  subroutine fill_field_with_data_real32_3d(fs, data, seed)
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: fs
    real(c_float), intent(inout) :: data(:,:,:)
    integer, intent(in) :: seed
    type(atlas_Field) :: field_global_index
    integer(ATLAS_KIND_GIDX), pointer :: global_index(:)
    integer(ATLAS_KIND_IDX) :: jblk, jpt, i
    integer :: j

    field_global_index = fs%global_index()
    call field_global_index%data(global_index)

    data = real(seed, c_float)
    do jblk = 1, fs%nblks()
      do i = 1, fs%block_size(jblk)
        jpt = fs%block_begin(jblk) + i - 1
        do j = 1, size(data,2)
          data(i,j,jblk) = real(seed, c_float) + real(global_index(jpt), c_float) + real(j, c_float)
        enddo
      enddo
    enddo
    call field_global_index%final()
  end subroutine fill_field_with_data_real32_3d

  subroutine fill_field_with_data_real64_4d(fs, data, seed)
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: fs
    real(c_double), intent(inout) :: data(:,:,:,:)
    integer, intent(in) :: seed
    type(atlas_Field) :: field_global_index
    integer(ATLAS_KIND_GIDX), pointer :: global_index(:)
    integer(ATLAS_KIND_IDX) :: jblk, jpt, i
    integer :: j, k

    field_global_index = fs%global_index()
    call field_global_index%data(global_index)

    data = real(seed, c_double)
    do jblk = 1, fs%nblks()
      do i = 1, fs%block_size(jblk)
        jpt = fs%block_begin(jblk) + i - 1
        do j = 1, size(data,2)
          do k = 1, size(data,3)
            data(i,j,k,jblk) = real(seed, c_double) + real(global_index(jpt), c_double) + real(j, c_double) + real(k, c_double)
          enddo
        enddo
      enddo
    enddo
    call field_global_index%final()
  end subroutine fill_field_with_data_real64_4d

  subroutine fill_field_with_data_real32_4d(fs, data, seed)
    type(atlas_functionspace_BlockStructuredColumns), intent(in) :: fs
    real(c_float), intent(inout) :: data(:,:,:,:)
    integer, intent(in) :: seed
    type(atlas_Field) :: field_global_index
    integer(ATLAS_KIND_GIDX), pointer :: global_index(:)
    integer(ATLAS_KIND_IDX) :: jblk, jpt, i
    integer :: j, k

    field_global_index = fs%global_index()
    call field_global_index%data(global_index)

    data = real(seed, c_float)
    do jblk = 1, fs%nblks()
      do i = 1, fs%block_size(jblk)
        jpt = fs%block_begin(jblk) + i - 1
        do j = 1, size(data,2)
          do k = 1, size(data,3)
            data(i,j,k,jblk) = real(seed, c_float) + real(global_index(jpt), c_float) + real(j, c_float) + real(k, c_float)
          enddo
        enddo
      enddo
    enddo
    call field_global_index%final()
  end subroutine fill_field_with_data_real32_4d

end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fcta_BlockStructuredColumns_checksum, fcta_BlockStructuredColumns_checksum_fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! ============================================================================
! Test: Checksum for single field (double precision, no levels)
! ============================================================================
TEST( test_blockstructuredcolumns_checksum_single_field_2d )
use atlas_functionspace_blockstructuredcolumns_module
use fckit_mpi_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field1, field2
character(len=256) :: checksum1, checksum2
type(fckit_mpi_comm) :: mpi

mpi = fckit_mpi_comm()

grid = atlas_StructuredGrid("O8")
fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=4)

! Create two identical fields and verify checksums are equal
field1 = fs%create_field(name="field1", kind=atlas_real(c_double))
field2 = fs%create_field(name="field2", kind=atlas_real(c_double))

call fill_field_with_data(fs, field1, 42)
call fill_field_with_data(fs, field2, 42)

checksum1 = fs%checksum(field1)
checksum2 = fs%checksum(field2)

if( mpi%rank() == 0 ) then
  write(msg,*) "Checksum 1: ", trim(checksum1)
  call atlas_log%info(msg)
  write(msg,*) "Checksum 2: ", trim(checksum2)
  call atlas_log%info(msg)
endif

FCTEST_CHECK_EQUAL( checksum1, checksum2 )

call field1%final()
call field2%final()
call fs%final()
call grid%final()
END_TEST

! ============================================================================
! Test: Checksum for single field (float precision, with levels)
! ============================================================================
TEST( test_blockstructuredcolumns_checksum_3d_field )
use atlas_functionspace_blockstructuredcolumns_module
use fckit_mpi_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field
character(len=256) :: checksum1, checksum2
type(fckit_mpi_comm) :: mpi

mpi = fckit_mpi_comm()

grid = atlas_StructuredGrid("O8")
fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=4, levels=5)

field = fs%create_field(name="field", kind=atlas_real(c_float))

call fill_field_with_data(fs, field, 123)

! First checksum
checksum1 = fs%checksum(field)

! Recompute should give same checksum
checksum2 = fs%checksum(field)

if( mpi%rank() == 0 ) then
  write(msg,*) "3D Field Checksum 1: ", trim(checksum1)
  call atlas_log%info(msg)
  write(msg,*) "3D Field Checksum 2: ", trim(checksum2)
  call atlas_log%info(msg)
endif

FCTEST_CHECK_EQUAL( checksum1, checksum2 )

call field%final()
call fs%final()
call grid%final()
END_TEST

! ============================================================================
! Test: Checksum for fieldset with multiple fields
! ============================================================================
TEST( test_blockstructuredcolumns_checksum_fieldset_multiple_fields )
use atlas_functionspace_blockstructuredcolumns_module
use fckit_mpi_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field1, field2, field3
type(atlas_FieldSet) :: fieldset, fieldset2
character(len=256) :: checksum1, checksum2, checksum_field1
type(fckit_mpi_comm) :: mpi

mpi = fckit_mpi_comm()

grid = atlas_StructuredGrid("O8")
fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=4)

! Create three fields with different data
field1 = fs%create_field(name="field1", kind=atlas_real(c_double))
field2 = fs%create_field(name="field2", kind=atlas_real(c_double))
field3 = fs%create_field(name="field3", kind=atlas_real(c_float))

call fill_field_with_data(fs, field1, 100)
call fill_field_with_data(fs, field2, 200)
call fill_field_with_data(fs, field3, 300)

! Create fieldsets with same fields in same order
fieldset = atlas_FieldSet()
call fieldset%add(field1)
call fieldset%add(field2)
call fieldset%add(field3)

fieldset2 = atlas_FieldSet()
call fieldset2%add(field1)
call fieldset2%add(field2)
call fieldset2%add(field3)

! Checksums of identical fieldsets should match
checksum1 = fs%checksum(fieldset)
checksum2 = fs%checksum(fieldset2)

if( mpi%rank() == 0 ) then
  write(msg,*) "FieldSet Checksum 1: ", trim(checksum1)
  call atlas_log%info(msg)
  write(msg,*) "FieldSet Checksum 2: ", trim(checksum2)
  call atlas_log%info(msg)
endif

FCTEST_CHECK_EQUAL( checksum1, checksum2 )

! Also verify that individual field checksum differs from combined fieldset
checksum_field1 = fs%checksum(field1)

if( mpi%rank() == 0 ) then
  write(msg,*) "Field 1 Checksum: ", trim(checksum_field1)
  call atlas_log%info(msg)
  write(msg,*) "FieldSet contains Field 1 but they have different checksums"
  call atlas_log%info(msg)
endif

! FieldSet checksum should differ from single field checksum (contains 3 fields)
FCTEST_CHECK( checksum1 /= checksum_field1 )

call fieldset%final()
call fieldset2%final()
call field1%final()
call field2%final()
call field3%final()
call fs%final()
call grid%final()
END_TEST

! ============================================================================
! Test: Checksum determinism with multiple variables per field
! ============================================================================
TEST( test_blockstructuredcolumns_checksum_multiple_variables )
use atlas_functionspace_blockstructuredcolumns_module
use fckit_mpi_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field1, field2
character(len=256) :: checksum1, checksum2
type(fckit_mpi_comm) :: mpi

mpi = fckit_mpi_comm()

grid = atlas_StructuredGrid("O8")
fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=4)

! Create fields with multiple variables
field1 = fs%create_field(name="multi_var1", kind=atlas_real(c_double), variables=3)
field2 = fs%create_field(name="multi_var2", kind=atlas_real(c_double), variables=3)

call fill_field_with_data(fs, field1, 555)
call fill_field_with_data(fs, field2, 555)

checksum1 = fs%checksum(field1)
checksum2 = fs%checksum(field2)

if( mpi%rank() == 0 ) then
  write(msg,*) "Multi-var Checksum 1: ", trim(checksum1)
  call atlas_log%info(msg)
  write(msg,*) "Multi-var Checksum 2: ", trim(checksum2)
  call atlas_log%info(msg)
endif

FCTEST_CHECK_EQUAL( checksum1, checksum2 )

call field1%final()
call field2%final()
call fs%final()
call grid%final()
END_TEST

! ============================================================================
! Test: Checksum invariance with halo - should exclude halo region
! Note: This test skipped because properly testing halo exclusion requires
! more careful setup of owned vs halo regions. The C++ tests verify this.
! ============================================================================
! TEST( test_blockstructuredcolumns_checksum_halo_consistency )
! SKIPPED - See C++ tests for halo invariance verification
!

! ============================================================================
! Test: Checksum differs when data changes
! ============================================================================
TEST( test_blockstructuredcolumns_checksum_data_sensitivity )
use atlas_functionspace_blockstructuredcolumns_module
use fckit_mpi_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field
character(len=256) :: checksum1, checksum2
real(c_double), pointer :: data_2d(:,:)
type(fckit_mpi_comm) :: mpi

mpi = fckit_mpi_comm()

grid = atlas_StructuredGrid("O8")
fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=4)

field = fs%create_field(name="field", kind=atlas_real(c_double))

! Fill with initial data
call fill_field_with_data(fs, field, 42)
checksum1 = fs%checksum(field)

! Modify data
call field%data(data_2d)
data_2d = data_2d + 1.0d0

checksum2 = fs%checksum(field)

if( mpi%rank() == 0 ) then
  write(msg,*) "Checksum before modification: ", trim(checksum1)
  call atlas_log%info(msg)
  write(msg,*) "Checksum after modification: ", trim(checksum2)
  call atlas_log%info(msg)
endif

! Checksums should differ after data modification
FCTEST_CHECK( checksum1 /= checksum2 )

call field%final()
call fs%final()
call grid%final()
END_TEST

! ============================================================================
! Test: Checksum determinism for different data types
! ============================================================================
TEST( test_blockstructuredcolumns_checksum_data_type_sensitivity )
use atlas_functionspace_blockstructuredcolumns_module
use fckit_mpi_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field_double, field_float, field1, field2
real(c_double), pointer :: data_dbl(:,:), data_dbl2(:,:)
real(c_float), pointer :: data_flt(:,:)
character(len=256) :: checksum_double, checksum_float, checksum1, checksum2
type(fckit_mpi_comm) :: mpi

mpi = fckit_mpi_comm()

grid = atlas_StructuredGrid("O16")
fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=8)

! ========================================================================
! Test 1: Fields created through functionspace, data wraps internal buffers
! (This demonstrates how BlockStructuredColumns manages field data internally)
! ========================================================================
field_double = fs%create_field(name="field_double", kind=atlas_real(c_double))
field_float = fs%create_field(name="field_float", kind=atlas_real(c_float))

! Access and fill the wrapped data
call field_double%data(data_dbl)
call field_float%data(data_flt)

! Fill with deterministic values
data_dbl = 1.618033988749895d0
data_flt = real(1.618033988749895d0, c_float)

! Compute checksums
checksum_double = fs%checksum(field_double)
checksum_float = fs%checksum(field_float)

if( mpi%rank() == 0 ) then
  write(msg,*) "Field wrapping double precision data: ", trim(checksum_double)
  call atlas_log%info(msg)
  write(msg,*) "Field wrapping single precision data: ", trim(checksum_float)
  call atlas_log%info(msg)
  write(msg,*) "Different data types produce different checksums"
  call atlas_log%info(msg)
endif

! Different data types should have different checksums
FCTEST_CHECK( checksum_double /= checksum_float )

! ========================================================================
! Test 2: Demonstrates that modifying the data changes checksum
! ========================================================================
field1 = fs%create_field(name="field1", kind=atlas_real(c_double))
field2 = fs%create_field(name="field2", kind=atlas_real(c_double))

call field1%data(data_dbl)
call field2%data(data_dbl2)

! Fill both with same data
data_dbl = 3.14159265358979d0
data_dbl2 = 3.14159265358979d0

checksum1 = fs%checksum(field1)
checksum2 = fs%checksum(field2)

if( mpi%rank() == 0 ) then
  write(msg,*) "Field 1 checksum (same data): ", trim(checksum1)
  call atlas_log%info(msg)
  write(msg,*) "Field 2 checksum (same data): ", trim(checksum2)
  call atlas_log%info(msg)
endif

FCTEST_CHECK_EQUAL( checksum1, checksum2 )

! Modify wrapped data in field1
data_dbl = 2.71828182845904d0

checksum1 = fs%checksum(field1)

if( mpi%rank() == 0 ) then
  write(msg,*) "Field 1 checksum (modified data): ", trim(checksum1)
  call atlas_log%info(msg)
  write(msg,*) "Checksum changed because wrapped data was modified"
  call atlas_log%info(msg)
endif

! Checksums should differ after modification
FCTEST_CHECK( checksum1 /= checksum2 )

call field_double%final()
call field_float%final()
call field1%final()
call field2%final()
call fs%final()
call grid%final()
END_TEST

! ============================================================================
! Test: Checksum consistency across rank types
! ============================================================================
TEST( test_blockstructuredcolumns_checksum_rank_consistency )
use atlas_functionspace_blockstructuredcolumns_module
use fckit_mpi_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field_rank2, field_rank3, field_rank4
character(len=256) :: checksum_r2, checksum_r3, checksum_r4
type(fckit_mpi_comm) :: mpi

mpi = fckit_mpi_comm()

grid = atlas_StructuredGrid("O8")
fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=4, levels=5)

! Create fields of different ranks
field_rank2 = fs%create_field(name="rank2_field", kind=atlas_real(c_double))
field_rank3 = fs%create_field(name="rank3_field", kind=atlas_real(c_double), levels=5)
field_rank4 = fs%create_field(name="rank4_field", kind=atlas_real(c_double), variables=2, levels=5)

call fill_field_with_data(fs, field_rank2, 111)
call fill_field_with_data(fs, field_rank3, 222)
call fill_field_with_data(fs, field_rank4, 333)

! Compute checksums for different rank fields
checksum_r2 = fs%checksum(field_rank2)
checksum_r3 = fs%checksum(field_rank3)
checksum_r4 = fs%checksum(field_rank4)

if( mpi%rank() == 0 ) then
  write(msg,*) "Rank 2 Field Checksum: ", trim(checksum_r2)
  call atlas_log%info(msg)
  write(msg,*) "Rank 3 Field Checksum: ", trim(checksum_r3)
  call atlas_log%info(msg)
  write(msg,*) "Rank 4 Field Checksum: ", trim(checksum_r4)
  call atlas_log%info(msg)
  write(msg,*) "All checksums should be different (different data and shapes)"
  call atlas_log%info(msg)
endif

! Checksums should be different (different data)
FCTEST_CHECK( checksum_r2 /= checksum_r3 )
FCTEST_CHECK( checksum_r3 /= checksum_r4 )
FCTEST_CHECK( checksum_r2 /= checksum_r4 )

call field_rank2%final()
call field_rank3%final()
call field_rank4%final()
call fs%final()
call grid%final()
END_TEST

TEST( test_blockstructuredcolumns_wrapped_fields )
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field_1, field_2, field_3, field_gfl
type(atlas_FieldSet) :: fieldset
real(c_double), allocatable :: gfl(:,:,:,:)
character(len=256) :: checksum
integer(c_int) :: nlev = 10, nproma = 8, nvar = 3

grid = atlas_StructuredGrid("O32")
fs = atlas_functionspace_BlockStructuredColumns(grid, nproma=nproma)
allocate(gfl(fs%nproma(),nlev,nvar,fs%nblks()))
call fill_field_with_data(fs, gfl, 42)

field_1 = atlas_Field(name="field_1", data=gfl(:,:,1,:))
checksum = fs%checksum(field_1)  ! Should compute checksum without error on wrapped data
write(msg,*) "Checksum for field_1 wrapping 3D slice: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "8b78")

field_2 = atlas_Field(name="field_2", data=gfl(:,:,2,:))
checksum = fs%checksum(field_2)  ! Should compute checksum without error on wrapped data
write(msg,*) "Checksum for field_2 wrapping 3D slice: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "31a0")

field_3 = atlas_Field(name="field_3", data=gfl(:,:,3,:))
checksum = fs%checksum(field_3)  ! Should compute checksum without error on wrapped data
write(msg,*) "Checksum for field_3 wrapping 3D slice: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "5496")

fieldset = atlas_FieldSet()
call fieldset%add(field_1)
call fieldset%add(field_2)
call fieldset%add(field_3)
checksum = fs%checksum(fieldset)  ! Should compute checksum without error on fieldset of wrapped fields
write(msg,*) "Checksum for fieldset of 3 wrapped fields: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "9305")

field_gfl = atlas_Field(name="gfl", data=gfl(:,:,:,:))
checksum = fs%checksum(field_gfl)  ! Should compute checksum without error on wrapped data
write(msg,*) "Checksum for field_gfl wrapping 4D array: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "338c")
END_TEST

TEST( test_blockstructuredcolumns_checksum_array_api )
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
real(c_double), allocatable :: gfl(:,:,:,:)
character(len=256) :: checksum
integer(c_int) :: nlev = 10, nproma = 8, nvar = 3

grid = atlas_StructuredGrid("O32")
fs = atlas_functionspace_BlockStructuredColumns(grid, nproma=nproma)
allocate(gfl(fs%nproma(),nlev,nvar,fs%nblks()))
call fill_field_with_data(fs, gfl, 42)

checksum = fs%checksum(gfl(:,:,1,:))
write(msg,*) "Checksum for array API wrapping 3D slice 1: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "8b78")

checksum = fs%checksum(gfl(:,:,2,:))
write(msg,*) "Checksum for array API wrapping 3D slice 2: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "31a0")

checksum = fs%checksum(gfl(:,:,3,:))
write(msg,*) "Checksum for array API wrapping 3D slice 3: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "5496")

checksum = fs%checksum(gfl(:,:,:,:))
write(msg,*) "Checksum for array API wrapping 4D array: ", trim(checksum)
call atlas_log%info(msg)
FCTEST_CHECK_EQUAL(checksum, "338c")

call fs%final()
call grid%final()
END_TEST

END_TESTSUITE
