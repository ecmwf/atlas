! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM TEST_SYNC_ON_FINAL
  USE ATLAS_MODULE, ONLY : ATLAS_LIBRARY
  IMPLICIT NONE

  CALL ATLAS_LIBRARY%INITIALISE()

  PRINT *, "begin 1: Without SYNC_ON_FINAL argument: default = .TRUE., MANUAL_SYNC not required"
  CALL RUN_TEST(MANUAL_SYNC=.TRUE.,  EXPECTED_FAILURE=.FALSE.)
  CALL RUN_TEST(MANUAL_SYNC=.FALSE., EXPECTED_FAILURE=.FALSE.)
  PRINT *, "end 1"
  PRINT *, ""

  PRINT *, "begin 2: With SYNC_ON_FINAL = .TRUE., MANUAL_SYNC not required"
  CALL RUN_TEST(SYNC_ON_FINAL=.TRUE., MANUAL_SYNC=.TRUE.,  EXPECTED_FAILURE=.FALSE.)
  CALL RUN_TEST(SYNC_ON_FINAL=.TRUE., MANUAL_SYNC=.FALSE., EXPECTED_FAILURE=.FALSE.)
  PRINT *, "end 2"
  PRINT *, ""

  PRINT *, "begin 3: With SYNC_ON_FINAL = .FALSE., MANUAL_SYNC required"
  CALL RUN_TEST(SYNC_ON_FINAL=.FALSE., MANUAL_SYNC=.TRUE.,  EXPECTED_FAILURE=.FALSE.)
  CALL RUN_TEST(SYNC_ON_FINAL=.FALSE., MANUAL_SYNC=.FALSE., EXPECTED_FAILURE=.TRUE.)
  PRINT *, "end 3"
  PRINT *, ""

  CALL ATLAS_LIBRARY%FINALISE()

CONTAINS

  SUBROUTINE RUN_TEST(SYNC_ON_FINAL, MANUAL_SYNC, EXPECTED_FAILURE)
    USE FIELD_MODULE
    USE PARKIND1
    USE FIELD_FACTORY_MODULE
    use iso_c_binding

    logical, intent(in), optional :: SYNC_ON_FINAL
    logical, intent(in) :: MANUAL_SYNC
    logical, intent(in) :: EXPECTED_FAILURE

    REAL(KIND=JPRB), ALLOCATABLE, TARGET :: D1(:,:,:,:,:), D2(:,:,:,:,:)
    CLASS(FIELD_4RB), POINTER :: W4 => NULL()
    REAL(KIND=JPRB), POINTER :: W4PTR(:,:,:,:)
    integer(kind=8) :: ptr

    ALLOCATE(D1(7, 9, 11, 13, 15))
    ALLOCATE(D2(7, 9, 11, 13, 15))
    D1 = 0
    D2 = 0

    if (.not. present(SYNC_ON_FINAL)) then
      CALL FIELD_NEW(W4, DATA=D1(1:1,:,:,:,3))
    else
      CALL FIELD_NEW(W4, DATA=D1(1:1,:,:,:,3), SYNC_ON_FINAL=SYNC_ON_FINAL)
    endif
    CALL W4%GET_HOST_DATA_RDWR(W4PTR)
    W4PTR=42
    CALL W4%GET_DEVICE_DATA_RDWR(W4PTR)
    !$ACC KERNELS DEFAULT(PRESENT)
    W4PTR(:,:,:,:)=92
    !$ACC END KERNELS
    D2(1:1,:,:,:,3)=92

    IF (MANUAL_SYNC) THEN
        CALL W4%GET_HOST_DATA_RDONLY(W4PTR)
    ENDIF

    CALL FIELD_DELETE(W4)

    IF (ANY(D1/=D2)) THEN
      IF (.NOT. EXPECTED_FAILURE) THEN
        write(0,*) "Unexpected failure: D1 /= D2"
        write(0,*) "D1 = ", D1(1:1,1,1,1,3)
        write(0,*) "D2 = ", D2(1:1,1,1,1,3)
        ERROR STOP
      ENDIF
    ENDIF
  END SUBROUTINE

END PROGRAM TEST_SYNC_ON_FINAL