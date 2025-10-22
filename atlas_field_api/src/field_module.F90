! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module field_module

#ifdef SINGLE_PRECISION

USE FIELD_RANKSUFF_MODULE, ONLY : &
  & FIELD_1RB => FIELD_1RM, &
  & FIELD_2RB => FIELD_2RM, &
  & FIELD_3RB => FIELD_3RM, &
  & FIELD_4RB => FIELD_4RM

#else

USE FIELD_RANKSUFF_MODULE, ONLY : &
  & FIELD_1RB => FIELD_1RD, &
  & FIELD_2RB => FIELD_2RD, &
  & FIELD_3RB => FIELD_3RD, &
  & FIELD_4RB => FIELD_4RD

#endif

USE FIELD_RANKSUFF_MODULE, ONLY : &
  & FIELD_1RD, FIELD_1RM, &
  & FIELD_2RD, FIELD_2RM, &
  & FIELD_3RD, FIELD_3RM, &
  & FIELD_4RD, FIELD_4RM

end module

