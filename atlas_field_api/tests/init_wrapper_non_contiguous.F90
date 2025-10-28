! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_WRAPPER_NON_CONTIGUOUS
        ! TEST IF WRAPPER WORKS WITH NON CONTIGUOUS ARRAY
        USE FIELD_ABORT_MODULE
        USE PARKIND1

        IMPLICIT NONE
        REAL(KIND=JPRB), ALLOCATABLE :: D(:,:)
        INTEGER :: I, J

        ALLOCATE(D(10,10))
        D=7
        CALL DO_STUFF_WITH_NONCONTIGUOUS_DATA(D(5:10,1:2))
        DO I=5,10
        DO J=1,2
        IF (D(I,J) /= 42) THEN
                CALL FIELD_ABORT ("ERROR")
        END IF
        END DO
        END DO

CONTAINS
        SUBROUTINE DO_STUFF_WITH_NONCONTIGUOUS_DATA(D)
                USE FIELD_MODULE
                USE FIELD_FACTORY_MODULE
                USE PARKIND1
                IMPLICIT NONE
                REAL(KIND=JPRB) :: D(:,:)
                CLASS(FIELD_2RB), POINTER :: W => NULL()
                REAL(KIND=JPRB), POINTER :: PTR(:,:)
                CALL FIELD_NEW(W, DATA=D)
                CALL W%GET_HOST_DATA_RDWR(PTR)
                PTR=42
                CALL FIELD_DELETE(W)
        END SUBROUTINE DO_STUFF_WITH_NONCONTIGUOUS_DATA

END PROGRAM INIT_WRAPPER_NON_CONTIGUOUS

