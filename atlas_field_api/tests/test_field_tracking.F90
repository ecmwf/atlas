! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM TEST_FIELD_TRACKING
  USE ATLAS_MODULE, ONLY : ATLAS_LIBRARY
  USE FIELD_MODULE
  USE FIELD_FACTORY_MODULE
  USE PARKIND1
  USE FIELD_ABORT_MODULE
  USE FIELD_UTIL_MODULE, ONLY : SET_DEVICE_MEMORY
  IMPLICIT NONE
  CLASS(FIELD_2RB), POINTER :: W => NULL()
  REAL(KIND=JPRB), ALLOCATABLE :: D(:,:)
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: D_GPU(:,:)
  REAL(KIND=JPRB), POINTER :: D_CPU(:,:)
  LOGICAL :: OKAY = .TRUE.
  INTEGER :: I,J

  CALL SET_DEVICE_MEMORY() ! Ensuring separate memory space for device, even in host-only mode

  CALL ATLAS_LIBRARY%INITIALISE()

  ALLOCATE(D(10,10))

  WRITE(0,*) "+ D = 42"
  D = 42._JPRB
  CALL FIELD_NEW(W, DATA=D)

  WRITE(0,*) "+ GET_DEVICE_DATA_WRONLY(D_GPU)  -- should not update device"
  CALL W%GET_DEVICE_DATA_WRONLY(D_GPU)

  WRITE(0,*) "+ ASSERT(D_GPU /= 42)"
  !$ACC SERIAL PRESENT (D_GPU) COPY(OKAY)
  DO J=1,10
    DO I=1,10
      IF(D_GPU(I,J) == 42._JPRB) THEN
        OKAY = .FALSE.
      END IF
    END DO
  END DO
  !$ACC END SERIAL
  IF (.NOT. OKAY) THEN
    CALL FIELD_ABORT("ERROR: Failed ASSERT(D_GPU /= 42)")
  END IF

  WRITE(0,*) "+ D_GPU = 84"
  !$ACC SERIAL PRESENT (D_GPU)
  DO J=1,10
    DO I=1,10
      D_GPU(I,J) = 84._JPRB
    END DO
  END DO
  !$ACC END SERIAL

  WRITE(0,*) "+ GET_HOST_DATA_RDONLY(D_CPU)    -- should update host"
  CALL W%GET_HOST_DATA_RDONLY(D_CPU)

  WRITE(0,*) "+ ASSERT(D_CPU == 84)"
  DO I=1,10
    DO J=1,10
      IF(D_CPU(I,J) /= 84._JPRB) THEN
        CALL FIELD_ABORT("ERROR: Failed ASSERT(D_CPU == 84)")
      END IF
    END DO
  END DO
  WRITE(0,*) "+ D_CPU = 21"
  D_CPU(1:10,1:10) = 21._JPRB

  WRITE(0,*) "+ GET_DEVICE_DATA_RDWR(D_GPU)    -- should not update device"
  CALL W%GET_DEVICE_DATA_RDWR(D_GPU)

  WRITE(0,*) "+ ASSERT(D_GPU == 84)"
  !$ACC SERIAL PRESENT (D_GPU) COPY(OKAY)
  DO I=1,10
    DO J=1,10
      IF(D_GPU(I,J) /= 84._JPRB) THEN
        OKAY = .FALSE.
      END IF
    END DO
  END DO
  !$ACC END SERIAL
  IF (.NOT. OKAY) THEN
    CALL FIELD_ABORT("ERROR: Failed ASSERT(D_GPU == 84)")
  END IF

  WRITE(0,*) "+ ASSERT(D_CPU == 21)"
  DO I=1,10
    DO J=1,10
      IF(D_CPU(I,J) /= 21._JPRB) THEN
          CALL FIELD_ABORT("ERROR: Failed ASSERT(D_CPU == 21)")
      END IF
    END DO
  END DO

  WRITE(0,*) "+ GET_HOST_DATA_RDONLY(D_CPU)    -- should update host"
  CALL W%GET_HOST_DATA_RDONLY(D_CPU)

  WRITE(0,*) "+ ASSERT(D_CPU == 84)"
  DO I=1,10
    DO J=1,10
      IF(D_CPU(I,J) /= 84._JPRB) THEN
          CALL FIELD_ABORT("ERROR: Failed ASSERT(D_CPU == 84)")
      END IF
    END DO
  END DO

  CALL FIELD_DELETE(W)

  CALL ATLAS_LIBRARY%FINALISE()

END PROGRAM TEST_FIELD_TRACKING

