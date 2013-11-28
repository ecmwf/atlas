! module split_globe_module
! Needs rtable.d file (can be found in ~rdx/data/ifs)

! This program is a stand-alone version of the algorithm used by IFS to split the Gaussian grid
! into regions for MPI parallelization. It is referred to as the "equal area" algorithm.

! The programs read the Gaussian Grid definition from file rtable.d.
! It then asks you how many regions you want to slit the globe in.
! The output is in file distribution.d.
! The output has one line per grid-point made up of 4 numbers:
! gridpoint number, latitude number, longitude number, and region
! The ordering of the grid-points is in the standard ECMWF GRIB order:
! from north to south and for each latitude from west to east starting at Greenwich

!=======================================================================

module eq_regions_mod
!
!     Purpose.
!     --------
!           eq_regions_mod provides the code to perform a high level 
!           partitioning of the surface of a sphere into regions of
!           equal area and small diameter.
!           the type.
!
!     Background.
!     -----------
!     This Fortran version of eq_regions is a much cut down version of the 
!     "Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox" of the
!     same name developed by Paul Leopardi at the University of New South Wales. 
!     This version has been coded specifically for the case of partitioning the 
!     surface of a sphere or S^dim (where dim=2) as denoted in the original code.
!     Only a subset of the original eq_regions package has been coded to determine 
!     the high level distribution of regions on a sphere, as the detailed 
!     distribution of grid points to each region is left to IFS software.
!     This is required to take into account the spatial distribution of grid 
!     points in an IFS gaussian grid and provide an optimal (i.e. exact) 
!     distribution of grid points over regions.
!
!     The following copyright notice for the eq_regions package is included from 
!     the original MatLab release.
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     + Release 1.10 2005-06-26                                                 +
!     +                                                                         +
!     + Copyright (c) 2004, 2005, University of New South Wales                 +
!     +                                                                         +
!     + Permission is hereby granted, free of charge, to any person obtaining   +
!     + a copy of this software and associated documentation files (the         +
!     + "Software"), to deal in the Software without restriction, including     +
!     + without limitation the rights to use, copy, modify, merge, publish,     +
!     + distribute, sublicense, and/or sell copies of the Software, and to      +
!     + permit persons to whom the Software is furnished to do so, subject to   +
!     + the following conditions:                                               +
!     +                                                                         +
!     + The above copyright notice and this permission notice shall be included +
!     + in all copies or substantial portions of the Software.                  +
!     +                                                                         +
!     + THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         +
!     + EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      +
!     + MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  +
!     + IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    +
!     + CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    +
!     + TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       +
!     + SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  +
!     +                                                                         +
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     Author.
!     -------
!        George Mozdzynski *ECMWF*
!
!     Modifications.
!     --------------
!        Original : 2006-04-15
!
!--------------------------------------------------------------------------------

use common_module  ,only : jpim     ,jprw

implicit none

save

private

public eq_regions,l_regions_debug,n_regions_ns,n_regions_ew,n_regions,my_region_ns,my_region_ew

real(kind=jprw) pi
logical :: l_regions_debug=.false.
integer(kind=jpim) :: n_regions_ns
integer(kind=jpim) :: n_regions_ew
integer(kind=jpim) :: my_region_ns
integer(kind=jpim) :: my_region_ew
integer(kind=jpim),allocatable :: n_regions(:)

CONTAINS

subroutine eq_regions(N)
!
! eq_regions uses the zonal equal area sphere partitioning algorithm to partition
! the surface of a sphere into N regions of equal area and small diameter.
!
integer(kind=jpim),intent(in) :: N
integer(kind=jpim) :: n_collars,j
real(kind=jprw),allocatable :: r_regions(:)
real(kind=jprw) :: c_polar

pi=2.0_jprw*asin(1.0_jprw)

n_regions(:)=0

if( N == 1 )then

  !
  ! We have only one region, which must be the whole sphere.
  !
  n_regions(1)=1
  n_regions_ns=1

else

  !
  ! Given N, determine c_polar
  ! the colatitude of the North polar spherical cap.
  !
  c_polar = polar_colat(N)
  !
  ! Given N, determine the ideal angle for spherical collars.
  ! Based on N, this ideal angle, and c_polar,
  ! determine n_collars, the number of collars between the polar caps.
  !
  n_collars = num_collars(N,c_polar,ideal_collar_angle(N))
  n_regions_ns=n_collars+2
  !
  ! Given N, c_polar and n_collars, determine r_regions,
  ! a list of the ideal real number of regions in each collar,
  ! plus the polar caps.
  ! The number of elements is n_collars+2.
  ! r_regions[1] is 1.
  ! r_regions[n_collars+2] is 1.
  ! The sum of r_regions is N.
  allocate(r_regions(n_collars+2))
  call ideal_region_list(N,c_polar,n_collars,r_regions)
  !
  ! Given N and r_regions, determine n_regions, a list of the natural number 
  ! of regions in each collar and the polar caps.
  ! This list is as close as possible to r_regions.
  ! The number of elements is n_collars+2.
  ! n_regions[1] is 1.
  ! n_regions[n_collars+2] is 1.
  ! The sum of n_regions is N.
  !
  call round_to_naturals(N,n_collars,r_regions)
  deallocate(r_regions)
  if( N /= sum(n_regions(:)) )then
    !write(*,'("eq_regions: N=",I10," sum(n_regions(:))=",I10)')N,sum(n_regions(:))
!    call abor1('eq_regions: N /= sum(n_regions)')
  endif

endif

if( l_regions_debug )then
  !write(*,'("eq_regions: N=",I6," n_regions_ns=",I4)') N,n_regions_ns
  do j=1,n_regions_ns
    !write(*,'("eq_regions: n_regions(",I4,")=",I4)') j,n_regions(j)
  enddo
endif
n_regions_ew=maxval(n_regions(:))

return
end subroutine eq_regions

function num_collars(N,c_polar,a_ideal) result(num_c)
!
!NUM_COLLARS The number of collars between the polar caps
!
! Given N, an ideal angle, and c_polar,
! determine n_collars, the number of collars between the polar caps.
!
integer(kind=jpim),intent(in) :: N
real(kind=jprw),intent(in) :: a_ideal,c_polar
integer(kind=jpim) :: num_c
logical enough
enough = (N > 2) .and. (a_ideal > 0)
if( enough )then
  num_c = max(1,nint((pi-2.*c_polar)/a_ideal))
else
  num_c = 0
endif
return
end function num_collars

subroutine ideal_region_list(N,c_polar,n_collars,r_regions)
!
!IDEAL_REGION_LIST The ideal real number of regions in each zone
!
! List the ideal real number of regions in each collar, plus the polar caps.
!
! Given N, c_polar and n_collars, determine r_regions, a list of the ideal real 
! number of regions in each collar, plus the polar caps.
! The number of elements is n_collars+2.
! r_regions[1] is 1.
! r_regions[n_collars+2] is 1.
! The sum of r_regions is N.
!
integer(kind=jpim),intent(in) :: N,n_collars
real(kind=jprw),intent(in) :: c_polar
real(kind=jprw),intent(out) :: r_regions(n_collars+2)
integer(kind=jpim) :: collar_n
real(kind=jprw) :: ideal_region_area,ideal_collar_area
real(kind=jprw) :: a_fitting
r_regions(:)=0.0_jprw
r_regions(1) = 1.0_jprw
if( n_collars > 0 )then
  !
  ! Based on n_collars and c_polar, determine a_fitting,
  ! the collar angle such that n_collars collars fit between the polar caps.
  !
  a_fitting = (pi-2.0_jprw*c_polar)/float(n_collars)
  ideal_region_area = area_of_ideal_region(N)
  do collar_n=1,n_collars
    ideal_collar_area = area_of_collar(c_polar+(collar_n-1)*a_fitting, &
     & c_polar+collar_n*a_fitting)
    r_regions(1+collar_n) = ideal_collar_area / ideal_region_area
  enddo
endif
r_regions(2+n_collars) = 1.
return
end subroutine ideal_region_list

function ideal_collar_angle(N) result(ideal)
!
! IDEAL_COLLAR_ANGLE The ideal angle for spherical collars of an EQ partition
!
! IDEAL_COLLAR_ANGLE(N) sets ANGLE to the ideal angle for the
! spherical collars of an EQ partition of the unit sphere S^2 into N regions.
!
integer(kind=jpim),intent(in) :: N
real(kind=jprw) :: ideal
ideal = area_of_ideal_region(N)**(0.5_jprw)
return
end function ideal_collar_angle

subroutine round_to_naturals(N,n_collars,r_regions)
!
! ROUND_TO_NATURALS Round off a given list of numbers of regions
!
! Given N and r_regions, determine n_regions, a list of the natural number 
! of regions in each collar and the polar caps.
! This list is as close as possible to r_regions, using rounding.
! The number of elements is n_collars+2.
! n_regions[1] is 1.
! n_regions[n_collars+2] is 1.
! The sum of n_regions is N.
!
integer(kind=jpim),intent(in) :: N,n_collars
real(kind=jprw),intent(in) :: r_regions(n_collars+2)
integer(kind=jpim) :: zone_n
real(kind=jprw) :: discrepancy
n_regions(1:n_collars+2) = r_regions(:)
discrepancy = 0.0_jprw
do zone_n = 1,n_collars+2
    n_regions(zone_n) = nint(r_regions(zone_n)+discrepancy);
    discrepancy = discrepancy+r_regions(zone_n)-float(n_regions(zone_n));
enddo
return
end subroutine round_to_naturals

function polar_colat(N) result(polar_c)
!
! Given N, determine the colatitude of the North polar spherical cap.
!
integer(kind=jpim),intent(in) :: N
real(kind=jprw) :: area
real(kind=jprw) :: polar_c
if( N == 1 ) polar_c=pi
if( N == 2 ) polar_c=pi/2.0_jprw
if( N > 2 )then
  area=area_of_ideal_region(N)
  polar_c=sradius_of_cap(area)
endif
return
end function polar_colat

function area_of_ideal_region(N) result(area)
!
! AREA_OF_IDEAL_REGION(N) sets AREA to be the area of one of N equal
! area regions on S^2, that is 1/N times AREA_OF_SPHERE.
!
integer(kind=jpim),intent(in) :: N
real(kind=jprw) :: area_of_sphere
real(kind=jprw) :: area
area_of_sphere = (2.0_jprw*pi**1.5_jprw/gamma(1.5_jprw))
area = area_of_sphere/float(N)
return
end function area_of_ideal_region

function sradius_of_cap(area) result(sradius)
!
! SRADIUS_OF_CAP(AREA) returns the spherical radius of
! an S^2 spherical cap of area AREA.
!
real(kind=jprw),intent(in) :: area
real(kind=jprw) :: sradius
sradius = 2.0_jprw*asin(sqrt(area/pi)/2.0_jprw)
return
end function sradius_of_cap

function area_of_collar(a_top, a_bot) result(area)
!
! AREA_OF_COLLAR Area of spherical collar
!
! AREA_OF_COLLAR(A_TOP, A_BOT) sets AREA to be the area of an S^2 spherical 
! collar specified by A_TOP, A_BOT, where A_TOP is top (smaller) spherical radius,
! A_BOT is bottom (larger) spherical radius.
!
real(kind=jprw),intent(in) :: a_top,a_bot
real(kind=jprw) area
area = area_of_cap(a_bot) - area_of_cap(a_top)
return
end function area_of_collar

function area_of_cap(s_cap) result(area)
!
! AREA_OF_CAP Area of spherical cap
!
! AREA_OF_CAP(S_CAP) sets AREA to be the area of an S^2 spherical
! cap of spherical radius S_CAP.
!
real(kind=jprw),intent(in) :: s_cap
real(kind=jprw) area
area = 4.0_jprw*pi * sin(s_cap/2.0_jprw)**2
return
end function area_of_cap

function gamma(x) result(gamma_res)
real(kind=jprw),intent(in) :: x
real(kind=jprw) :: gamma_res
real(kind=jprw) :: p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
real(kind=jprw) :: w,y
integer(kind=jpim) :: k,n
parameter (&
& p0 =   0.999999999999999990e+00_jprw,&
& p1 =  -0.422784335098466784e+00_jprw,&
& p2 =  -0.233093736421782878e+00_jprw,&
& p3 =   0.191091101387638410e+00_jprw,&
& p4 =  -0.024552490005641278e+00_jprw,&
& p5 =  -0.017645244547851414e+00_jprw,&
& p6 =   0.008023273027855346e+00_jprw)
parameter (&
& p7 =  -0.000804329819255744e+00_jprw,&
& p8 =  -0.000360837876648255e+00_jprw,&
& p9 =   0.000145596568617526e+00_jprw,&
& p10 = -0.000017545539395205e+00_jprw,&
& p11 = -0.000002591225267689e+00_jprw,&
& p12 =  0.000001337767384067e+00_jprw,&
& p13 = -0.000000199542863674e+00_jprw)
n = nint(x - 2)
w = x - (n + 2)
y = ((((((((((((p13 * w + p12) * w + p11) * w + p10) *&
&    w + p9) * w + p8) * w + p7) * w + p6) * w + p5) *&
&    w + p4) * w + p3) * w + p2) * w + p1) * w + p0
if (n .gt. 0) then
  w = x - 1
  do k = 2, n
    w = w * (x - k)
  end do
else
  w = 1
  do k = 0, -n - 1
    y = y * (x + k)
  end do
end if
gamma_res = w / y
return
end function gamma

end module eq_regions_mod

!=======================================================================

module SUSTAONL_MOD
contains
SUBROUTINE SUSTAONL(KMEDIAP,KRESTM,LDWEIGHTED_DISTR,PWEIGHT,PMEDIAP,KPROCAGP, &
 & NPROC,KDGL,KDLON,KLOEN,LDSPLIT,KFRSTLAT,KLSTLAT,KPTRFRSTLAT,KMYSETA,KSTA,KONL)

!**** *SUSTAONL * - Routine to initialize parallel environment

!     Purpose.
!     --------
!           Initialize NSTA and NONL.
!           Calculation of distribution of grid points to processors :
!           Splitting of grid in B direction

!**   Interface.
!     ----------
!        *CALL* *SUSTAONL *

!        Explicit arguments : 
!        --------------------
!                     KMEDIAP    - mean number of grid points per PE
!                     KRESTM     - number of PEs with one extra point
!                     LDWEIGHTED_DISTR -true if weighted distribution
!                     PWEIGHT    -weight per grid-point if weighted distribution
!                     PMEDIAP    -mean weight per PE if weighted distribution
!                     KPROCAGP   -number of grid points per A set

!        Implicit arguments :
!        --------------------


!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!        Modified 98-08-10 by K. YESSAD: removal of LRPOLE option.
!          - removal of LRPOLE in YOMCT0.
!          - removal of code under LRPOLE.
!        Modified 98-12-04 C. Fischer: merge with SUESTAONL (Aladin)
!        R. El Khatib 05-Apr-2007 Enable back vectorization on NEC
!     ------------------------------------------------------------------

USE common_module  ,ONLY : JPIM     ,jprw
USE EQ_REGIONS_MOD

IMPLICIT NONE


!     DUMMY 
INTEGER(KIND=JPIM),INTENT(IN) :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(IN) :: KRESTM
REAL(KIND=jprw),INTENT(IN)    :: PWEIGHT(:)
LOGICAL,INTENT(IN)            :: LDWEIGHTED_DISTR
REAL(KIND=jprw),INTENT(IN)    :: PMEDIAP
INTEGER(KIND=JPIM),INTENT(IN) :: KPROCAGP(:)
INTEGER(KIND=JPIM),INTENT(IN) :: NPROC
INTEGER(KIND=JPIM),INTENT(IN) :: KDGL
INTEGER(KIND=JPIM),INTENT(IN) :: KDLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLOEN(:)
LOGICAL           ,INTENT(IN) :: LDSPLIT
INTEGER(KIND=JPIM),INTENT(IN) :: KFRSTLAT(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KLSTLAT(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KPTRFRSTLAT(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KMYSETA
INTEGER(KIND=JPIM),INTENT(OUT) :: KSTA(:,:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KONL(:,:)
!     LOCAL

INTEGER(KIND=JPIM) :: IXPTLAT(KDGL), ILSTPTLAT(KDGL),ISENDREQ(NPROC)
INTEGER(KIND=JPIM) :: ICHK(KDLON,KDGL), ICOMBUF(KDGL*N_REGIONS_EW*2)
INTEGER(KIND=JPIM) :: I1, I2, IBUFLEN, IDGLG, IDWIDE,&
             &IGL, IGL1, IGL2, IGLOFF, IGPTA, &
             &IGPTPRSETS, IGPTS, IGPTSP, ILEN, ILRECV, &
             &ILSEND, INPLAT, INXLAT, IPOS, &
             &IPROCB, IPTSRE, IRECV, IPE, &
             &IREST, ISEND, ITAG, JA, JB, JGL, JL, JNPTSRE, &
             &ILAT, ILON, ILOEN
INTEGER(KIND=JPIM),ALLOCATABLE :: ICOMBUFG(:)
REAL(KIND=jprw),ALLOCATABLE :: ZWEIGHT(:,:)
INTEGER(KIND=JPIM) :: JJ, ILENG(NPROC), IOFF(NPROC)

LOGICAL :: LLABORT
LOGICAL :: LLP1,LLP2

REAL(KIND=jprw) ::  ZLAT, ZLAT1, ZCOMP
REAL(KIND=jprw) :: ZDIVID(KDGL),ZXPTLAT(KDGL)

!      -----------------------------------------------------------------

IXPTLAT  (:)=999999
ILSTPTLAT(:)=999999

LLP1 = .TRUE.
LLP2 = .TRUE.

IDWIDE  = KDGL/2
IBUFLEN = KDGL*N_REGIONS_EW*2
IDGLG   = KDGL

I1 = MAX(   1,KFRSTLAT(KMYSETA)-(KFRSTLAT(KMYSETA)-1))
I2 = MIN(IDGLG,KLSTLAT (KMYSETA)-(KFRSTLAT(KMYSETA)-1))

ILEN = KLSTLAT(KMYSETA) - KFRSTLAT(KMYSETA)+1

IGPTPRSETS = SUM(KLOEN(1:KFRSTLAT(KMYSETA)-1))

IF (LDSPLIT) THEN
  IGPTA=0
  DO JA=1,KMYSETA-1
    IGPTA = IGPTA + KPROCAGP(JA)
  ENDDO
  IGPTS = KPROCAGP(KMYSETA)
ELSE
  IGPTA = IGPTPRSETS
  IGPTS = SUM(KLOEN(KFRSTLAT(KMYSETA):KLSTLAT(KMYSETA)))
ENDIF

IGPTSP = IGPTS/N_REGIONS(KMYSETA)
IREST = IGPTS-N_REGIONS(KMYSETA)*IGPTSP
IXPTLAT(1) = IGPTA-IGPTPRSETS+1
ZXPTLAT(1) = REAL(IXPTLAT(1))
ILSTPTLAT(1) = KLOEN(KFRSTLAT(KMYSETA))
INPLAT = KLOEN(KFRSTLAT(KMYSETA))-IXPTLAT(1)+1
DO JGL=2,ILEN
  IXPTLAT(JGL) = 1
  ZXPTLAT(JGL) = 1.0_jprw
  ILSTPTLAT(JGL) =  KLOEN(KFRSTLAT(KMYSETA)+JGL-1)
  INPLAT = INPLAT+KLOEN(KFRSTLAT(KMYSETA)+JGL-1)
ENDDO
ILSTPTLAT(ILEN) = KLOEN(KLSTLAT(KMYSETA))-INPLAT+IGPTS


!  grid point decomposition
!  ---------------------------------------
DO JGL=1,ILEN
  ZDIVID(JGL)=REAL(KLOEN(KFRSTLAT(KMYSETA)+JGL-1),jprw)
ENDDO
IF( LDWEIGHTED_DISTR )THEN
  ALLOCATE(ZWEIGHT(KLOEN(KDGL/2),KDGL))
  IGL=0
  DO JGL=1,KDGL
    DO JL=1,KLOEN(JGL)
      IGL=IGL+1
      ZWEIGHT(JL,JGL)=PWEIGHT(IGL)
    ENDDO
  ENDDO
  ZCOMP=0
  IGPTS=0
ENDIF

DO JB=1,N_REGIONS(KMYSETA)

  IF( .NOT.LDWEIGHTED_DISTR )THEN

    IF (JB <= IREST) THEN
      IPTSRE = IGPTSP+1
    ELSE
      IPTSRE = IGPTSP
    ENDIF

    DO JNPTSRE=1,IPTSRE

      ZLAT  = 1._jprw
      ZLAT1 = 1._jprw

      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          ZLAT1 = (ZXPTLAT(JGL)-1.0_jprw)/ZDIVID(JGL)
          IF (ZLAT1 < ZLAT) THEN
            ZLAT   = ZLAT1
            INXLAT = JGL
          ENDIF
        ENDIF
      ENDDO
  
      IF (INXLAT >= I1 .AND. INXLAT <= I2) THEN
        IGL=(KPTRFRSTLAT(KMYSETA)-1)+INXLAT
        IF (KSTA(IGL,JB) == 0) THEN
          KSTA(IGL,JB) = IXPTLAT(INXLAT)
        ENDIF
        KONL(IGL,JB) = KONL(IGL,JB)+1
      ENDIF
      IXPTLAT(INXLAT) = IXPTLAT(INXLAT)+1
      ZXPTLAT(INXLAT) = REAL(IXPTLAT(INXLAT),jprw)
    ENDDO

  ELSE

    DO WHILE ( (JB <  N_REGIONS(KMYSETA) .AND. ZCOMP < PMEDIAP) &
        & .OR. (JB == N_REGIONS(KMYSETA) .AND. IGPTS < KPROCAGP(KMYSETA)) )

      IGPTS = IGPTS + 1
      ZLAT  = 1._jprw
      ZLAT1 = 1._jprw

      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          ZLAT1 = (ZXPTLAT(JGL)-1.0_jprw)/ZDIVID(JGL)
          IF (ZLAT1 < ZLAT) THEN
            ZLAT   = ZLAT1
            INXLAT = JGL
          ENDIF
        ENDIF
      ENDDO
  
      IF (INXLAT >= I1 .AND. INXLAT <= I2) THEN
        IGL=(KPTRFRSTLAT(KMYSETA)-1)+INXLAT
        IF (KSTA(IGL,JB) == 0) THEN
          KSTA(IGL,JB) = IXPTLAT(INXLAT)
        ENDIF
        KONL(IGL,JB) = KONL(IGL,JB)+1
        IF(IGL<1.OR.IGL>KDGL+N_REGIONS_NS-1)THEN
          !write(0,*) 'SUSTAONL: IGL<1.OR.IGL>KDGL+N_REGIONS_NS-1'
          stop 1
        ENDIF
        ILON=KSTA(IGL,JB)+KONL(IGL,JB)-1
        ILAT=KFRSTLAT(KMYSETA)+INXLAT-1
        ILOEN=KLOEN(ILAT)
        IF(ILON<1.OR.ILON>ILOEN)THEN
          !write(0,*)' SUSTAONL: ILON<1.OR.ILON>ILOEN'
          stop 1
        ENDIF
        ZCOMP = ZCOMP + ZWEIGHT(ILON,ILAT)
      ENDIF
      IXPTLAT(INXLAT) = IXPTLAT(INXLAT)+1
      ZXPTLAT(INXLAT) = REAL(IXPTLAT(INXLAT),jprw)
    ENDDO

    ZCOMP = ZCOMP - PMEDIAP

  ENDIF

ENDDO

IF( LDWEIGHTED_DISTR )THEN
  DEALLOCATE(ZWEIGHT)
ENDIF


IF(KMYSETA == N_REGIONS_NS) THEN
  LLABORT = .FALSE.
  DO JGL=1,KDGL
    DO JL=1,KLOEN(JGL)
      ICHK(JL,JGL) = 1
    ENDDO
  ENDDO
  DO JA=1,N_REGIONS_NS
    IGLOFF = KPTRFRSTLAT(JA)
    DO JB=1,N_REGIONS(JA)
      IGL1 = KFRSTLAT(JA)
      IGL2 = KLSTLAT(JA)
      DO JGL=IGL1,IGL2
        IGL = IGLOFF+JGL-IGL1
        DO JL=KSTA(IGL,JB),KSTA(IGL,JB)+KONL(IGL,JB)-1
          IF( ICHK(JL,JGL) /= 1 )THEN
            !WRITE(0,'(" SUSTAONL : seta=",i4," setb=",i4,&
            ! &" row=",I4," sta=",I4," INVALID GRID POINT")')&
            ! &JA,JB,JGL,JL
            !WRITE(0,'(" SUSTAONL : seta=",i4," setb=",i4,&
            ! &" ROW=",I4," sta=",I4," INVALID GRID POINT")')&
            ! &JA,JB,JGL,JL
            LLABORT = .TRUE.
          ENDIF
          ICHK(JL,JGL) = 2
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO JGL=1,KDGL
    DO JL=1,KLOEN(JGL)
      IF( ICHK(JL,JGL) /= 2 )THEN
        !WRITE(0,'(" SUSTAONL : row=",i4," sta=",i4,&
        ! &" GRID POINT NOT ASSIGNED")') JGL,JL
        LLABORT = .TRUE.
      ENDIF
    ENDDO
  ENDDO
  IF( LLABORT )THEN
    !WRITE(0,'(" SUSTAONL : inconsistent partitioning")')
    STOP 1
  ENDIF

  IF (LLP1) THEN
    !WRITE(UNIT=6,FMT='('' OUTPUT FROM ROUTINE SUSTAONL '')')
    !WRITE(UNIT=6,FMT='('' '')')
    !WRITE(UNIT=6,FMT='('' PARTITIONING INFORMATION '')')
    !WRITE(UNIT=6,FMT='('' '')')
    IPROCB = MIN(32,N_REGIONS_EW)
    !WRITE(UNIT=6,FMT='(17X," SETB=",32(1X,I3))') (JB,JB=1,IPROCB)
    DO JA=1,N_REGIONS_NS
      IPROCB = MIN(32,N_REGIONS(JA))
      !WRITE(UNIT=6,FMT='('' '')')
      IGLOFF = KPTRFRSTLAT(JA)
      IGL1 = KFRSTLAT(JA)
      IGL2 = KLSTLAT(JA)
      DO JGL=IGL1,IGL2
        IGL=IGLOFF+JGL-IGL1
        !WRITE(UNIT=6,FMT='(" SETA=",I3," LAT=",I3," NSTA=",&
        ! &32(1X,I3))') JA,JGL,(KSTA(IGL,JB),JB=1,IPROCB)
        !WRITE(UNIT=6,FMT='(" SETA=",I3," LAT=",I3," KONL=",&
        ! &32(1X,I3))') JA,JGL,(KONL(IGL,JB),JB=1,IPROCB)
        !WRITE(UNIT=6,FMT='('' '')')
      ENDDO
      !WRITE(UNIT=6,FMT='('' '')')
    ENDDO
    !WRITE(UNIT=6,FMT='('' '')')
    !WRITE(UNIT=6,FMT='('' '')')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE SUSTAONL

END MODULE SUSTAONL_MOD

!=======================================================================

module SUMPLATBEQ_MOD
contains
SUBROUTINE SUMPLATBEQ(KDGL,KPROC,KPROCA,KLOENG,LDSPLIT,&
                    &PWEIGHT,LDWEIGHTED_DISTR,PMEDIAP,KPROCAGP,&
                    &KMEDIAP,KRESTM,KINDIC,KLAST)

!**** *SUMPLATBEQ * - Routine to initialize parallel environment
!                     (latitude partitioning for LEQ_REGIONS=T)

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *SUMPLATBEQ *

!     Explicit arguments - input :
!     -------------------- 
!                          KDGL       -last  latitude
!                          KPROC      -total number of processors
!                          KPROCA     -number of processors in A direction
!                          KLOENG     -actual number of longitudes per latitude.
!                          LDSPLIT    -true for latitudes shared between sets
!                          PWEIGHT    -weight per grid-point if weighted distribution
!                          LDWEIGHTED_DISTR -true if weighted distribution

!     Explicit arguments - output:
!     -------------------- 
!                          PMEDIAP    -mean weight per PE if weighted distribution
!                          KMEDIAP    -mean number of grid points per PE
!                          KPROCAGP   -number of grid points per A set
!                          KRESTM     -number of PEs with one extra point
!                          KINDIC     -intermediate quantity for 'sumplat'
!                          KLAST      -intermediate quantity for 'sumplat'

!        Implicit arguments :
!        --------------------


!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        G. Mozdzynski

!     Modifications.
!     --------------
!        Original : April 2006
!     ------------------------------------------------------------------


USE common_module  ,ONLY : JPIM     ,jprw

USE EQ_REGIONS_MOD

IMPLICIT NONE


!     * DUMMY:
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROC
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCA
INTEGER(KIND=JPIM),INTENT(IN)  :: KLOENG(KDGL)
REAL(KIND=jprw),INTENT(IN)     :: PWEIGHT(:)
LOGICAL,INTENT(IN)  :: LDSPLIT
LOGICAL,INTENT(INOUT)  :: LDWEIGHTED_DISTR
REAL(KIND=jprw),INTENT(OUT)     :: PMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT)  :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT)  :: KRESTM
INTEGER(KIND=JPIM),INTENT(OUT)  :: KINDIC(KPROCA)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KLAST(KPROCA)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KPROCAGP(KPROCA)

!     * LOCAL:

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ICOMP, IGL, IMAXI, IMEDIA, IMEDIAP, ITOT, JA, JB, IA, JGL,&
            &ILAST,IREST,IPE,I2REGIONS,IGP
REAL(KIND=jprw) :: ZMEDIA, ZCOMP
LOGICAL   :: LLDONE

!      -----------------------------------------------------------------

!*       1.    COMPUTATION OF KMEDIAP, KRESTM, KINDIC, KLAST.
!              ----------------------------------------------
100 CONTINUE
!     * Computation of KMEDIAP and KRESTM.

IF (.NOT.LDWEIGHTED_DISTR) THEN

  IMEDIA = SUM(KLOENG(1:KDGL))
  KMEDIAP = IMEDIA / KPROC

  IF( KPROC > 1 )THEN
! test if KMEDIAP is too small and no more than 2 asets would be required 
! for the first latitude
    IF( LDSPLIT )THEN
      I2REGIONS=N_REGIONS(1)+N_REGIONS(2)
      IF( KMEDIAP < (KLOENG(1)-1)/I2REGIONS+1 )THEN
        !WRITE(0,'("SUMPLATBEQ: KMEDIAP=",I6," I2REGIONS=",I3," KLOENG(1)=",I4)')&
        !&KMEDIAP,I2REGIONS,KLOENG(1)
!        CALL ABORT_TRANS ('SUMPLATBEQ: NPROC TOO BIG FOR THIS RESOLUTION, LDSPLIT=T')
      ENDIF
    ELSE
! test for number asets too large for the number of latitudes
      IF( KPROCA > KDGL )THEN
        !WRITE(0,'("SUMPLATBEQ: KMEDIAP=",I6," KPROCA=",I4," KDGL=",I4)')&
        !&KMEDIAP,KPROCA,KDGL
 !       CALL ABORT_TRANS ('SUMPLATBEQ: NPROC TOO BIG FOR THIS RESOLUTION, LDSPLIT=F')
      ENDIF
    ENDIF
  ENDIF

  KRESTM = IMEDIA - KMEDIAP * KPROC
  IF (KRESTM  >  0) KMEDIAP = KMEDIAP + 1

ELSE

  ZMEDIA = SUM(PWEIGHT(:))
  PMEDIAP = ZMEDIA / KPROC

ENDIF

!     * Computation of intermediate quantities KINDIC and KLAST

IF (LDSPLIT) THEN

  KPROCAGP(:)=0
  IREST = 0
  ILAST =0
  IPE=0
  ZCOMP=0
  IGP=0
  DO JA=1,KPROCA
    ICOMP=0
    DO JB=1,N_REGIONS(JA)
      IF( LDWEIGHTED_DISTR )THEN
        DO WHILE ( ( JA == KPROCA .OR. ZCOMP < PMEDIAP ) .AND. IGP < SIZE(PWEIGHT) )
          IGP = IGP + 1
          ICOMP = ICOMP + 1
          ZCOMP = ZCOMP + PWEIGHT(IGP)
        ENDDO
        ZCOMP = ZCOMP - PMEDIAP
      ELSE
        IPE=IPE+1
        IF (IPE  <=  KRESTM .OR. KRESTM  ==  0) THEN
          ICOMP = ICOMP + KMEDIAP
        ELSE
          ICOMP = ICOMP + (KMEDIAP-1)
        ENDIF
      ENDIF
    ENDDO
    KPROCAGP(JA)=ICOMP
    ITOT = IREST
    IGL = ILAST+1
    DO JGL=IGL,KDGL
      ILAST = JGL
      IF(ITOT+KLOENG(JGL) < ICOMP) THEN
        ITOT = ITOT+KLOENG(JGL)
      ELSEIF(ITOT+KLOENG(JGL) == ICOMP) THEN
        IREST = 0
        KLAST(JA) = JGL 
        KINDIC(JA) = 0
        EXIT
      ELSE
        IREST =  KLOENG(JGL) -(ICOMP-ITOT)
        KLAST(JA) = JGL 
        KINDIC(JA) = JGL
        EXIT
      ENDIF
    ENDDO
  ENDDO
  IF( LDWEIGHTED_DISTR )THEN
    IF( KLAST(KPROCA) /= KDGL )THEN
      DO JA=1,KPROCA
!        IF( MYPROC == 1 )THEN
          !WRITE(0,'("SUMPLATBEQ_MOD: JA=",I3," KLAST=",I3," KINDIC=",I3)')&
          !&JA,KLAST(JA),KINDIC(JA)
!        ENDIF
      ENDDO
      !WRITE(0,'("SUMPLATBEQ: LWEIGHTED_DISTR=T FAILED TO PARTITION GRID, REVERTING TO ",&
      !& " LWEIGHTED_DISTR=F PARTITIONING")')
      LDWEIGHTED_DISTR=.FALSE.
      GOTO 100
    ENDIF
  ENDIF
  IF( SUM(KPROCAGP(:)) /= SUM(KLOENG(1:KDGL)) )THEN
    !WRITE(0,'("SUM(KPROCAGP(:))=",I12)')SUM(KPROCAGP(:))
    !WRITE(0,'("SUM(KLOENG(:))=",I12)')SUM(KLOENG(1:KDGL))
    STOP 1
  ENDIF

ELSE

  IF( LDWEIGHTED_DISTR )THEN
    !WRITE(0,*)'SUMPLATBEQ: LSPLIT=F NOT SUPPORTED FOR WEIGHTED DISTRIBUTION '
    STOP 1
  ENDIF

  KINDIC(:) = 0
  LLDONE=.FALSE.
  IMEDIAP=KMEDIAP
  !WRITE(0,'("SUMPLATBEQ: IMEDIAP=",I6)')IMEDIAP
  DO WHILE(.NOT.LLDONE)
!   loop until a satisfactory distribution can be found
    IA=1
    IMAXI=IMEDIAP*N_REGIONS(IA)
    DO JGL=1,KDGL
      KLAST(IA)=JGL
      IMAXI=IMAXI-KLOENG(JGL)
      IF( IA == KPROCA .AND. JGL == KDGL )THEN
        !WRITE(0,'("SUMPLATBEQ: EXIT 1")')
        EXIT
      ENDIF
      IF( IA == KPROCA .AND. JGL < KDGL )THEN
        !WRITE(0,'("SUMPLATBEQ: EXIT 2")')
        KLAST(KPROCA)=KDGL
        EXIT
      ENDIF
      IF( IA < KPROCA .AND. JGL == KDGL )THEN
        DO JA=KPROCA,IA+1,-1
          KLAST(JA)=KDGL+JA-KPROCA
        ENDDO
        DO JA=KPROCA,2,-1
          IF( KLAST(JA) <= KLAST(JA-1) )THEN
            KLAST(JA-1)=KLAST(JA)-1
          ENDIF
        ENDDO
         ! WRITE(0,'("SUMPLATBEQ: EXIT 3")')
        EXIT
      ENDIF
      IF( IMAXI <= 0 )THEN
        IA=IA+1
        IMAXI=IMAXI+IMEDIAP*N_REGIONS(IA)
      ENDIF
    ENDDO
    IF( KPROCA > 1 .AND. KLAST(KPROCA) == KLAST(KPROCA-1) )THEN
      IMEDIAP=IMEDIAP-1
        !WRITE(0,'("SUMPLATBEQ: REDUCING IMEDIAP=",I6)')IMEDIAP
      IF( IMEDIAP <= 0 )THEN
        !WRITE(0,*) 'SUMPLATBEQ: PROBLEM PARTITIONING WITH LSPLIT=F, IMEDIAP <= 0'
        STOP 1
      ENDIF
    ELSE
      LLDONE=.TRUE.
    ENDIF
  ENDDO
    
ENDIF

END SUBROUTINE SUMPLATBEQ

END MODULE SUMPLATBEQ_MOD


!===========================================================================

module split_globe_module
  use common_module, only: jprw, jpim, log_str, log_info
  implicit none
  private

  public split_globe

contains

  subroutine split_globe(rtable,nproc,proc,glb_idx)
    character(len=*), intent(in)    :: rtable   ! rtable.d file
    integer(kind=jpim), intent(in)  :: nproc    ! Number of partitions
    integer(kind=jpim), intent(out) :: proc(:)
    integer(kind=jpim), intent(out) :: glb_idx(:)
    integer(kind=jpim)              :: ndgl     ! Number of lattitudes
    integer(kind=jpim), allocatable :: nloen(:) ! Number of longitude points for each latitude
    integer(kind=jpim) :: jlat,idum
    open(21,file=rtable,access='sequential',status='unknown')
    read(21,*)
    read(21,*)ndgl
    allocate(nloen(ndgl))
    do jlat=1,ndgl
      read(21,*) idum,nloen(jlat)
    enddo
    close(21)
    !write(log_str,*) 'number of latitudes,total number of grid-points ',ndgl,sum(nloen); call log_info()
    call split_points(nproc,ndgl,nloen,proc,glb_idx)
  end subroutine split_globe


  SUBROUTINE SPLIT_POINTS(NPROC,NDGL,NLOEN,proc,glb_idx)
    USE EQ_REGIONS_MOD
    USE SUMPLATBEQ_MOD
    USE SUSTAONL_MOD
    IMPLICIT NONE
    INTEGER(KIND=JPIM),INTENT(IN) :: NPROC        ! Number of tasks (regions)
    INTEGER(KIND=JPIM),INTENT(IN) :: NDGL         ! Number of latitude rows
    INTEGER(KIND=JPIM),INTENT(IN) :: NLOEN(NDGL)  ! Number of longitude points for each latitude
    INTEGER(KIND=JPIM),INTENT(OUT):: proc(:)      ! Return processor ranks for each point
    INTEGER(KIND=JPIM),INTENT(OUT):: glb_idx(:)      ! Return processor ranks for each point


    INTEGER(KIND=JPIM) JNS,JEW,JLAT,IDUM
    INTEGER(KIND=JPIM),ALLOCATABLE :: INDIC(:),ILAST(:),IPROCAGP(:)
    INTEGER(KIND=JPIM),ALLOCATABLE :: IFRSTLAT(:),ILSTLAT(:),IPTRLAT(:)
    INTEGER(KIND=JPIM),ALLOCATABLE :: IPTRFRSTLAT(:),IPTRLSTLAT(:)
    INTEGER(KIND=JPIM),ALLOCATABLE :: IONL(:,:),ISTA(:,:),IGLOBAL(:,:)
    INTEGER(KIND=JPIM),ALLOCATABLE :: NGLOBALPROC(:)
    INTEGER(KIND=JPIM) :: IPTRLATITUDE,IDLON,IND,IPROC,IGL,I,JL
    INTEGER(KIND=JPIM) :: IMEDIAP,IRESTM
    REAL(KIND=jprw) :: ZMEDIAP,ZWEIGHT(1)
    LOGICAL :: LLSPLIT,LLWEIGHTED_DISTR
    LOGICAL,ALLOCATABLE :: LLSPLITLAT(:)

    ALLOCATE(N_REGIONS(NPROC+2))
    CALL EQ_REGIONS(NPROC)

    !print *,N_REGIONS_NS, 'NORTH-SOUTH BANDS'
    !print *,'EAST-WEST BANDS'
    !DO JNS=1,N_REGIONS_NS
    !  PRINT *,JNS,N_REGIONS(JNS)
    !ENDDO

    LLSPLIT=.TRUE.
    LLWEIGHTED_DISTR=.FALSE.
    ALLOCATE(INDIC(N_REGIONS_NS))
    ALLOCATE(ILAST(N_REGIONS_NS))
    ALLOCATE(IPROCAGP(N_REGIONS_NS))
    CALL SUMPLATBEQ(NDGL,NPROC,N_REGIONS_NS,NLOEN,LLSPLIT,ZWEIGHT,LLWEIGHTED_DISTR , &
     & ZMEDIAP,IPROCAGP,IMEDIAP,IRESTM,INDIC,ILAST)

    ALLOCATE(IFRSTLAT(N_REGIONS_NS))
    ALLOCATE(ILSTLAT(N_REGIONS_NS))
    IFRSTLAT(1) = 1
    ILSTLAT(N_REGIONS_NS) = NDGL
    DO JNS=1,N_REGIONS_NS-1
      IF ((.NOT. LLSPLIT) .OR. INDIC(JNS) == 0) THEN
        IFRSTLAT(JNS+1) = ILAST(JNS) + 1
        ILSTLAT(JNS) = ILAST(JNS)
      ELSE
        IFRSTLAT(JNS+1) = INDIC(JNS)
        ILSTLAT(JNS) = INDIC(JNS)
      ENDIF
    ENDDO
    !PRINT *,'IFRSTLAT ',IFRSTLAT
    !PRINT *,'ILSTLAT ',ILSTLAT

    ALLOCATE(IPTRLAT(NDGL))
    ALLOCATE(LLSPLITLAT(NDGL))
    DO JLAT=1,NDGL
      IPTRLAT  (JLAT)=-999
      LLSPLITLAT(JLAT)=.FALSE.
    ENDDO
    IPTRLATITUDE=0
    DO JNS=1,N_REGIONS_NS
      DO JLAT=IFRSTLAT(JNS),ILSTLAT(JNS)
        IPTRLATITUDE=IPTRLATITUDE+1
        LLSPLITLAT(JLAT)=.TRUE.
        IF( IPTRLAT(JLAT) == -999 )THEN
          IPTRLAT(JLAT)=IPTRLATITUDE
          LLSPLITLAT(JLAT)=.FALSE.
        ENDIF
      ENDDO
    ENDDO
    ALLOCATE(IPTRFRSTLAT(N_REGIONS_NS))
    ALLOCATE(IPTRLSTLAT(N_REGIONS_NS))
    DO JNS=1,N_REGIONS_NS
      IF( LLSPLITLAT(IFRSTLAT(JNS)) .AND. JNS /= 1)THEN
        IPTRFRSTLAT(JNS)=IPTRLAT(IFRSTLAT(JNS))+1
      ELSE
        IPTRFRSTLAT(JNS)=IPTRLAT(IFRSTLAT(JNS))
      ENDIF
      IF( LLSPLITLAT(ILSTLAT(JNS)) .AND. JNS == N_REGIONS_NS)THEN
        IPTRLSTLAT(JNS)=IPTRLAT(ILSTLAT(JNS))+1
      ELSE
        IPTRLSTLAT(JNS)=IPTRLAT(ILSTLAT(JNS))
      ENDIF
    ENDDO

    ALLOCATE(ISTA(NDGL+N_REGIONS_NS-1,N_REGIONS_EW))
    ALLOCATE(IONL(NDGL+N_REGIONS_NS-1,N_REGIONS_EW))
    ISTA=0
    IONL=0
    IDLON=MAXVAL(NLOEN)
    DO JNS=1,N_REGIONS_NS
        CALL SUSTAONL(IMEDIAP,IRESTM,LLWEIGHTED_DISTR,ZWEIGHT,ZMEDIAP,IPROCAGP, &
         & NPROC,NDGL,IDLON,NLOEN,LLSPLIT,IFRSTLAT,ILSTLAT,IPTRFRSTLAT,JNS,ISTA,IONL)
    ENDDO


    ALLOCATE(IGLOBAL(IDLON,NDGL))
    I=0
    IGLOBAL(:,:)=0
    DO JLAT=1,NDGL
      DO JL=1,NLOEN(JLAT)
        I=I+1
        IGLOBAL(JL,JLAT) = I
      ENDDO
    ENDDO

    ALLOCATE(NGLOBALPROC(SUM(NLOEN(:))))
    IPROC=0
    DO JNS=1,N_REGIONS_NS
      DO JEW=1,N_REGIONS(JNS)
        IPROC=IPROC+1
        I=0
        DO JLAT=IFRSTLAT(JNS),ILSTLAT(JNS)
          IGL = IPTRFRSTLAT(JNS)+JLAT-IFRSTLAT(JNS)
          DO JL=ISTA(IGL,JEW),ISTA(IGL,JEW)+IONL(IGL,JEW)-1
            IND=IGLOBAL(JL,JLAT)
            NGLOBALPROC(IND)=IPROC
            I=I+1
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    I=0
    DO JLAT=1,NDGL

      DO JL=1,NLOEN(JLAT)
        I=I+1
        IND=IGLOBAL(JL,JLAT)
        proc(I) = NGLOBALPROC(IND)-1
        glb_idx(I) = IND
      ENDDO

      ! Periodic point duplicated for every lattitude
      JL = 1
      I=I+1
      IND=IGLOBAL(JL,JLAT)
      glb_idx(I) = IND
      proc(I) = -1 ! invalid, marked for removal NGLOBALPROC(IND)-1
    ENDDO
    CLOSE(22)

  END SUBROUTINE SPLIT_POINTS

end module split_globe_module
