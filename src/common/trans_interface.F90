module trans_interface
use common_module, only : jprb
use MPL_module, only: MPL_INIT, MPL_NUMPROC
implicit none
save
private

public trans_setup
public trans_filter
public trans_zero_divergence

integer, parameter :: output_unit = 0
logical :: lfirst_call = .True.

#include "setup_trans0.h"
#include "setup_trans.h"
#include "dir_trans.h"
#include "inv_trans.h"
#include "dist_grid.h"
#include "gath_spec.h"
#include "trans_inq.h"
#include "specnorm.h"

contains

subroutine trans_setup(rtable,IGPTOT,handle)
  character(len=*), intent(in) :: rtable
  integer, intent(out) :: IGPTOT, handle

  integer              :: ndgl     ! Number of lattitudes
  integer, allocatable :: nloen(:) ! Number of longitude points for each latitude
  integer :: jlat,idum,NSMAX,nproc

  open(21,file=rtable,access='sequential',status='unknown')
  read(21,*)
  read(21,*)ndgl
  allocate(nloen(ndgl))
  do jlat=1,ndgl
    read(21,*) idum,nloen(jlat)
  enddo
  close(21)

  NSMAX = NDGL-1

  if (lfirst_call) then
    CALL MPL_INIT()
    nproc = MPL_NUMPROC
    CALL SETUP_TRANS0(KOUT=output_unit,KERR=0,KPRINTLEV=0,KMAX_RESOL=10,&
      &               KPRGPNS=NPROC,KPRGPEW=1,KPRTRW=NPROC,LDEQ_REGIONS=.True.)
    lfirst_call = .False.
  end if
  CALL SETUP_TRANS(KSMAX=NSMAX,KDGL=NDGL,KLOEN=NLOEN,LDSPLIT=.True.,&
    &              KRESOL=handle)

  CALL TRANS_INQ(KRESOL=handle,KGPTOT=IGPTOT)

end subroutine trans_setup

subroutine trans_norm(handle,field2d,norm)
  integer, intent(in) :: handle
  real(kind=jprb), intent(in) :: field2d(:,:,:)
  real(kind=jprb), allocatable, intent(out) :: norm(:)
  real(kind=jprb), allocatable :: ZSPEC(:,:)
  integer :: IMAXFLD, NSPEC2, IGPTOT
  CALL TRANS_INQ(KRESOL=handle,KSPEC2=NSPEC2)
  CALL TRANS_INQ(KRESOL=handle,KGPTOT=IGPTOT)
  IMAXFLD = size(field2d,2)
  ALLOCATE(ZSPEC(IMAXFLD,NSPEC2))
  if( .not. allocated(norm) ) ALLOCATE(NORM(IMAXFLD))

  CALL DIR_TRANS(KRESOL=handle,PSPSCALAR=ZSPEC,PGP=field2d,KPROMA=IGPTOT)

  CALL SPECNORM(KRESOL=handle,PSPEC=ZSPEC,PNORM=norm)
end subroutine trans_norm

subroutine trans_filter(handle,field2d,cutoff)
  integer, intent(in) :: handle
  real(kind=jprb), intent(inout) :: field2d(:,:,:)
  real(kind=jprb), intent(in) :: cutoff
  real(kind=jprb), allocatable :: ZSPEC(:,:)
  integer :: NPROMA, IGPTOT, NSPEC2,NSPEC2G,IMAXFLD,NUMP,NSMAX,IFILT
  integer :: IFILT, JMLOC,IM,JN,IJSE,JFLD
  integer, allocatable :: MYMS(:), NASM0(:)
  CALL TRANS_INQ(KRESOL=handle,KGPTOT=IGPTOT)
  CALL TRANS_INQ(KRESOL=handle,KSPEC2=NSPEC2,KSPEC2G=NSPEC2G)
  IMAXFLD = size(field2d,2)
  ALLOCATE(ZSPEC(IMAXFLD,NSPEC2))
  NPROMA=IGPTOT

  !field2d (:,:,:) = 1.

  CALL DIR_TRANS(KRESOL=handle,PSPSCALAR=ZSPEC,PGP=field2d,KPROMA=NPROMA)

  if (cutoff /= 1._jprb) then
    CALL TRANS_INQ(KNUMP=NUMP,KSMAX=NSMAX)
    allocate( myms(nump) )
    CALL TRANS_INQ(KMYMS=MYMS)

    allocate( NASM0(0:NSMAX) )
    CALL TRANS_INQ(KASM0=NASM0)

    !NUMP,MYMS,NASM0 from TRANS_INQ
    IFILT=NSMAX*cutoff+1
    DO JMLOC=1,NUMP
      IM=MYMS(JMLOC)
      DO JN=MAX(IM,IFILT),NSMAX
        DO JFLD=1,IMAXFLD
          IJSE=NASM0(IM)+2*(JN-IM)
          ZSPEC(JFLD,IJSE  ) = 0.0
          ZSPEC(JFLD,IJSE+1) = 0.0
        ENDDO
      ENDDO
    ENDDO
  end if
  CALL INV_TRANS(KRESOL=handle,PSPSCALAR=ZSPEC,PGP=field2d,KPROMA=NPROMA)
end subroutine trans_filter

subroutine trans_zero_divergence(handle,field2d)
  integer, intent(in) :: handle
  real(kind=jprb), intent(inout) :: field2d(:,:,:)
  real(kind=jprb), allocatable :: ZVOR(:,:), ZDIV(:,:), ZSPEC(:,:)
  integer :: NPROMA, IGPTOT, NSPEC2,NSPEC2G,IMAXFLD
  CALL TRANS_INQ(KRESOL=handle,KGPTOT=IGPTOT)
  CALL TRANS_INQ(KRESOL=handle,KSPEC2=NSPEC2,KSPEC2G=NSPEC2G)
  IMAXFLD = size(field2d,2)
  ALLOCATE(ZVOR(IMAXFLD/2,NSPEC2))
  ALLOCATE(ZDIV(IMAXFLD/2,NSPEC2))
  ALLOCATE(ZSPEC(IMAXFLD,NSPEC2))

  NPROMA=IGPTOT

  if (size(field2d,1) /= IGPTOT) then
    write(0,*) "nb_nodes is different: ", shape(field2d), IGPTOT
  end if

  !CALL DIR_TRANS(KRESOL=handle,PSPVOR=ZVOR,PSPDIV=ZDIV,PGP=field2d,KPROMA=NPROMA)
  CALL DIR_TRANS(KRESOL=handle,PSPSCALAR=ZSPEC,PGP=field2d,KPROMA=NPROMA)

  !ZDIV(:,:) = 0.

  !CALL INV_TRANS(KRESOL=handle,PSPVOR=ZVOR,PSPDIV=ZDIV,PGP=field2d,KPROMA=NPROMA)
  CALL INV_TRANS(KRESOL=handle,PSPSCALAR=ZSPEC,PGP=field2d,KPROMA=NPROMA)

end subroutine trans_zero_divergence

end module trans_interface