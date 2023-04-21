!----------------------------------------------------------------------
! Module: stripack
!> STRIPACK routines
! Source: https://people.sc.fsu.edu/~jburkardt/f_src/stripack/stripack.html
! Author: Robert Renka
! Original licensing: none
! Modified by Benjamin Menetrier for BUMP
! Modified by Francois Hebert for OOPS
! Licensing: this code is distributed under the CeCILL-C license
! Copyright 2022- UCAR
!----------------------------------------------------------------------

module stripack_mod

implicit none

integer, parameter :: kind_real = 8

real(kind_real), parameter :: zero    = 0.0_kind_real
real(kind_real), parameter :: one     = 1.0_kind_real
real(kind_real), parameter :: pi      = acos(-1.0_kind_real)
real(kind_real), parameter :: four    = 4.0_kind_real
real(kind_real), parameter :: hundred = 100.0_kind_real

interface addnod
   module procedure stripack_addnod
end interface
interface area
   module procedure stripack_area
end interface
interface bdyadd
   module procedure stripack_bdyadd
end interface
interface bnodes
   module procedure stripack_bnodes
end interface
interface covsph
   module procedure stripack_covsph
end interface
interface det
   module procedure stripack_det
end interface
interface insert
   module procedure stripack_insert
end interface
interface intadd
   module procedure stripack_intadd
end interface
interface jrand
   module procedure stripack_jrand
end interface
interface left
   module procedure stripack_left
end interface
interface lstptr
   module procedure stripack_lstptr
end interface
interface swap
   module procedure stripack_swap
end interface
interface swptst
   module procedure stripack_swptst
end interface
interface trfind
   module procedure stripack_trfind
end interface
interface trmesh
   module procedure stripack_trmesh
end interface
interface trlist2
   module procedure stripack_trlist2
end interface

private
public :: area,bnodes,trfind,trmesh,trlist2

contains

subroutine abor1_ftn( msg )
implicit none
character(len=*), intent(in) :: msg
!write(0,*) msg
end subroutine

!----------------------------------------------------------------------
! For quick testing, FH copied here a few dependencies from SABER repo
! TODO implement a cleaner solution to avoid duplication across repos
!----------------------------------------------------------------------

function indist(x,y) result(test)
implicit none
real(kind_real), intent(in) :: x !< First real
real(kind_real), intent(in) :: y !< Second real
real(kind_real), parameter ::rth = 1.0e-12
logical :: test
test = .true.
if ((abs(x)>zero).or.(abs(y)>zero)) then
  test = abs(x-y)<rth*(abs(x+y))
end if
end function indist

function eq(x,y) result(test)
implicit none
real(kind_real), intent(in) :: x !< First real
real(kind_real), intent(in) :: y !< Second real
logical :: test
test = indist(x,y)
end function eq

function inf(x,y) result(test)
implicit none
real(kind_real), intent(in) :: x !< First real
real(kind_real), intent(in) :: y !< Second real
logical :: test
test = (x<y)
test = test.and.(.not.indist(x,y))
end function inf

function sup(x,y) result(test)
implicit none
real(kind_real), intent(in) :: x !< First real
real(kind_real), intent(in) :: y !< Second real
logical :: test
test = (x>y)
test = test.and.(.not.indist(x,y))
end function sup

function supeq(x,y) result(test)
implicit none
real(kind_real), intent(in) :: x !< First real
real(kind_real), intent(in) :: y !< Second real
logical :: test
test = sup(x,y).or.eq(x,y)
end function supeq

function small(x,y) result(test)
implicit none
real(kind_real), intent(in) :: x !< First real
real(kind_real), intent(in) :: y !< Second real
real(kind_real), parameter ::rth = 1.0e-12
logical :: test
test = abs(x)<rth*abs(y)
end function small


!----------------------------------------------------------------------
! Subroutine: stripack_addnod
!> Add a node to a triangulation
!----------------------------------------------------------------------
subroutine stripack_addnod(nst,k,x,y,z,list,lptr,lend,lnew,ier)

! Passed variables
integer,intent(in) :: nst           !< Index of a node at which TRFIND, begins its search. Search
                                    !< time depends on the proximity of this node to K. If NST<1,
                                    !< the search is begun at node K-1.
integer,intent(in) :: k             !< Nodal index (index for X, Y, Z, and LEND) of the new node to
                                    !< be added. 4<=K.
real(kind_real),intent(in) :: x(k)  !< X-coordinates of the nodes
real(kind_real),intent(in) :: y(k)  !< Y-coordinates of the nodes
real(kind_real),intent(in) :: z(k)  !< Z-coordinates of the nodes
integer,intent(inout) :: list(*)    !< On input, the data structure associated with the
                                    !< triangulation of nodes 1 to K-1. On output, the data has been
                                    !< updated to include node K. The array lengths are assumed to
                                    !< be large enough to add node K. Refer to TRMESH.
integer,intent(inout) :: lptr(*)    !< cf. list
integer,intent(inout) :: lend(k)    !< cf. list
integer,intent(inout) :: lnew       !< cf. list
integer,intent(out)   :: ier        !< Error return code
                                    !<     0 if no errors were encountered.
                                    !<    -1 if K is outside its valid range on input.
                                    !<    -2 if all nodes (including K) are collinear (lie on a common geodesic).
                                    !<     L if nodes L and K coincide for some L < K.


! Local variables
integer :: i1,i2,i3,in1,io1,io2,ist,kk,km1,l,lp,lpf,lpo1,lpo1s
real(kind_real) :: b1,b2,b3,p(3)
logical :: output

! Initialization
ier = 0
kk = k
if (kk<4) then
   call abor1_ftn('stripack_addnod: K<4')
   ier = -1
   return
endif
km1 = kk-1
ist = nst
if (ist<1) ist = km1
p(1) = x(kk)
p(2) = y(kk)
p(3) = z(kk)

! Find a triangle (I1,I2,I3) containing K or the rightmost (I1) and leftmost (I2) visible boundary
! nodes as viewed from node K.
call trfind(ist,p,km1,x,y,z,list,lptr,lend,b1,b2,b3,i1,i2,i3)

! Test for colinear or duplicate nodes.
if (i1==0) then
   call abor1_ftn('stripack_addnod: colinear points')
   ier = -2
   return
endif
if (i3/=0) then
  l = i1
  if (eq(p(1),x(l)).and.eq(p(2),y(l)).and.eq(p(3),z(l))) then
    call abor1_ftn('stripack_addnod: colinear or duplicate nodes')
    ier = l
    return
  end if
  l = i2
  if (eq(p(1),x(l)).and.eq(p(2),y(l)).and.eq(p(3),z(l))) then
    call abor1_ftn('stripack_addnod: colinear or duplicate nodes')
    ier = l
    return
  end if
  l = i3
  if (eq(p(1),x(l)).and.eq(p(2),y(l)).and.eq(p(3),z(l))) then
    call abor1_ftn('stripack_addnod: colinear or duplicate nodes')
    ier = l
    return
  end if
  call intadd(kk,i1,i2,i3,list,lptr,lend,lnew)
else
  if (i1/=i2) then
    call bdyadd(kk,i1,i2,list,lptr,lend,lnew)
  else
    call covsph(kk,i1,list,lptr,lend,lnew)
  end if
end if

! Initialize variables for optimization of the triangulation.
lp = lend(kk)
lpf = lptr(lp)
io2 = list(lpf)
lpo1 = lptr(lpf)
io1 = abs(list(lpo1))

! Begin loop: find the node opposite K.
do
  call lstptr(lend(io1),io2,list,lptr,lp)
  if (0<=list(lp)) then
    lp = lptr(lp)
    in1 = abs(list(lp))

    ! Swap test:  if a swap occurs,two new arcs are opposite K and must be tested.
    lpo1s = lpo1
    call swptst(in1,kk,io1,io2,x,y,z,output)
    if (.not.output) then
      if ((lpo1==lpf).or.(list(lpo1)<0)) exit
      io2 = io1
      lpo1 = lptr(lpo1)
      io1 = abs(list(lpo1))
      cycle
    end if
    call swap(in1,kk,io1,io2,list,lptr,lend,lpo1)

    ! A swap is not possible because KK and IN1 are already adjacent. This error in SWPTST only
    ! occurs in the neutral case and when there are nearly duplicate nodes.
    if (lpo1/=0) then
      io1 = in1
      cycle
    end if
    lpo1 = lpo1s
  end if

  ! No swap occurred. Test for termination and reset IO2 and IO1.
  if ((lpo1==lpf).or.(list(lpo1)<0)) exit
  io2 = io1
  lpo1 = lptr(lpo1)
  io1 = abs(list(lpo1))
end do

end subroutine stripack_addnod

!----------------------------------------------------------------------
! Subroutine: stripack_area
!> Compute the area of a spherical triangle
!----------------------------------------------------------------------
subroutine stripack_area(v1,v2,v3,area)

! Passed variables
real(kind_real),intent(in) :: v1(3) !< First point cartesian coordinates
real(kind_real),intent(in) :: v2(3) !< Second point cartesian coordinates
real(kind_real),intent(in) :: v3(3) !< Third point cartesian coordinates
real(kind_real),intent(out) :: area !< Area on the unit sphere

! Local variables
real(kind_real) :: a1,a2,a3,ca1,ca2,ca3,s12,s23,s31
real(kind_real) :: dv1(3),dv2(3),dv3(3),u12(3),u23(3),u31(3)

! Initialization
dv1 = v1
dv2 = v2
dv3 = v3

! Compute cross products Uij = Vi X Vj
u12(1) = dv1(2)*dv2(3)-dv1(3)*dv2(2)
u12(2) = dv1(3)*dv2(1)-dv1(1)*dv2(3)
u12(3) = dv1(1)*dv2(2)-dv1(2)*dv2(1)
u23(1) = dv2(2)*dv3(3)-dv2(3)*dv3(2)
u23(2) = dv2(3)*dv3(1)-dv2(1)*dv3(3)
u23(3) = dv2(1)*dv3(2)-dv2(2)*dv3(1)
u31(1) = dv3(2)*dv1(3)-dv3(3)*dv1(2)
u31(2) = dv3(3)*dv1(1)-dv3(1)*dv1(3)
u31(3) = dv3(1)*dv1(2)-dv3(2)*dv1(1)

! Normalize Uij to unit vectors
s12 = dot_product(u12,u12)
s23 = dot_product(u23,u23)
s31 = dot_product(u31,u31)

! Test for a degenerate triangle associated with collinear vertices
if ((.not.(abs(s12)>zero)).or.(.not.(abs(s23)>zero)).or.(.not.(abs(s31)>zero))) then
   area = zero
else
   s12 = sqrt(s12)
   s23 = sqrt(s23)
   s31 = sqrt(s31)
   u12 = u12/s12
   u23 = u23/s23
   u31 = u31/s31

   ! Compute interior angles Ai as the dihedral angles between planes
   ca1 = -dot_product(u12,u31)
   ca2 = -dot_product(u23,u12)
   ca3 = -dot_product(u31,u23)
   ca1 = max(ca1,-one)
   ca1 = min(ca1,+one)
   ca2 = max(ca2,-one)
   ca2 = min(ca2,+one)
   ca3 = max(ca3,-one)
   ca3 = min(ca3,+one)
   a1 = acos(ca1)
   a2 = acos(ca2)
   a3 = acos(ca3)

   ! Compute area = a1 + a2 + a3 - pi
   area = a1+a2+a3-pi
   if (area<zero) area = zero
end if

end subroutine stripack_area


!----------------------------------------------------------------------
! Subroutine: stripack_bdyadd
!> Add a boundary node to a triangulation
!----------------------------------------------------------------------
subroutine stripack_bdyadd(kk,i1,i2,list,lptr,lend,lnew)

! Passed variable
integer,intent(in) :: kk         !< Index of a node to be connected to the sequence of all visible
                                 !< boundary nodes. 1<=KK and KK must not be equal to I1 or I2.
integer,intent(in) :: i1         !< First (rightmost as viewed from KK) boundary node in the
                                 !< triangulation that is visible from node KK (the line segment
                                 !< KK-I1 intersects no arcs).
integer,intent(in) :: i2         !< Last (leftmost) boundary node that is visible from node KK. I1
                                 !< and I2 may be determined by TRFIND.
integer,intent(inout) :: list(*) !< Triangulation data structure created by TRMESH. Nodes I1 and I2
                                 !< must be included in the triangulation. On output, the data
                                 !< structure is updated with the addition of node KK. Node KK is
                                 !< connected to I1, I2, and all boundary nodes in between.
integer,intent(inout) :: lptr(*) !< cf. list
integer,intent(inout) :: lend(*) !< cf. list
integer,intent(inout) :: lnew    !< cf. list

! Local variables
integer :: k,lp,lsav,n1,n2,next,nsav

! Initialization
k = kk
n1 = i1
n2 = i2

! Add K as the last neighbor of N1.
lp = lend(n1)
lsav = lptr(lp)
lptr(lp) = lnew
list(lnew) = -k
lptr(lnew) = lsav
lend(n1) = lnew
lnew = lnew+1
next = -list(lp)
list(lp) = next
nsav = next

! Loop on the remaining boundary nodes between N1 and N2, adding K as the first neighbor.
do
  lp = lend(next)
  call insert(k,lp,list,lptr,lnew)
  if (next==n2) exit
  next = -list(lp)
  list(lp) = next
end do

! Add the boundary nodes between N1 and N2 as neighbors of node K.
lsav = lnew
list(lnew) = n1
lptr(lnew) = lnew+1
lnew = lnew+1
next = nsav
do
  if (next==n2) exit
  list(lnew) = next
  lptr(lnew) = lnew+1
  lnew = lnew+1
  lp = lend(next)
  next = list(lp)
end do
list(lnew) = -n2
lptr(lnew) = lsav
lend(k) = lnew
lnew = lnew+1

end subroutine stripack_bdyadd

!----------------------------------------------------------------------
! Subroutine: stripack_bnodes
!> Return the boundary nodes of a triangulation
!----------------------------------------------------------------------
subroutine stripack_bnodes(n,list,lptr,lend,nodes,nb)

! Passed variables
integer,intent(in) :: n             !< Number of nodes in the triangulation. 3 <= N.
integer,intent(in) :: list(6*(n-2)) !< The data structure defining the triangulation, created by
                                    !< TRMESH.
integer,intent(in) :: lptr(6*(n-2)) !< cf. list
integer,intent(in) :: lend(n)       !< cf. list
integer,intent(out) :: nodes(n)     !< Ordered sequence of NB boundary node indexes in the range 1
                                    !< to N.
integer,intent(out) :: nb           !< Number of boundary nodes.

! Local variables
integer :: i,k,lp,n0,na,nn,nst,nt

! Initialization
nn = n

! Search for a boundary node
nst = 0
do i = 1, nn
   lp = lend(i)
   if (list(lp)<0) then
      nst = i
      exit
   end if
end do

! The triangulation contains no boundary nodes
if (nst==0) then
   nb = 0
   na = 3*(nn-2)
   nt = 2*(nn-2)
   return
end if

! NST is the first boundary node encountered
! Initialize for traversal of the boundary
nodes(1) = nst
k = 1
n0 = nst

! Traverse the boundary in counterclockwise order
do
   lp = lend(n0)
   lp = lptr(lp)
   n0 = list(lp)
   if (n0==nst) then
      exit
   end if
   k = k+1
   nodes(k) = n0
end do

! Store the counts
nb = k
nt = 2*n-nb-2
na = nt+n-1

end subroutine stripack_bnodes

!----------------------------------------------------------------------
! Subroutine: stripack_covsph
!> Connect an exterior node to boundary nodes, covering the sphere
!----------------------------------------------------------------------
subroutine stripack_covsph(kk,n0,list,lptr,lend,lnew)

! Passed variables
integer,intent(in) :: kk         !< Index of the node to be connected to the set of all boundary
                                 !< nodes. 4<=KK.
integer,intent(in) :: n0         !< Index of a boundary node (in the range 1 to KK-1). N0 may be
                                 !< determined by TRFIND.
integer,intent(inout) :: list(*) !< Triangulation data structure created by TRMESH. Node N0 must be
                                 !< included in the triangulation. On output, updated with the
                                 !< addition of node KK as the last entry. The updated triangulation
                                 !< contains no boundary nodes.
integer,intent(inout) :: lptr(*) !< cf. list
integer,intent(inout) :: lend(*) !< cf. list
integer,intent(inout) :: lnew    !< cf. list

! Local variables
integer :: k,lp,lsav,next,nst

! Initialization
k = kk
nst = n0

! Traverse the boundary in clockwise order, inserting K as the first neighbor of each boundary node,
! and converting the boundary node to an interior node.
next = nst
do
   lp = lend(next)
   call insert(k,lp,list,lptr,lnew)
   next = -list(lp)
   list(lp) = next
   if (next==nst) exit
end do

! Traverse the boundary again, adding each node to K's adjacency list.
lsav = lnew
do
   lp = lend(next)
   list(lnew) = next
   lptr(lnew) = lnew+1
   lnew = lnew+1
   next = list(lp)
  if (next==nst) exit
end do
lptr(lnew-1) = lsav
lend(k) = lnew-1

end subroutine stripack_covsph

!----------------------------------------------------------------------
! Subroutine: stripack_det
!> Compute 3D determinant
!----------------------------------------------------------------------
subroutine stripack_det(x1,y1,z1,x2,y2,z2,x0,y0,z0,output)

! Passed variables
real(kind_real),intent(in) :: x1      !< X-coordinate, term 1
real(kind_real),intent(in) :: y1      !< Y-coordinate, term 1
real(kind_real),intent(in) :: z1      !< Z-coordinate, term 1
real(kind_real),intent(in) :: x2      !< X-coordinate, term 2
real(kind_real),intent(in) :: y2      !< Y-coordinate, term 2
real(kind_real),intent(in) :: z2      !< Z-coordinate, term 2
real(kind_real),intent(in) :: x0      !< X-coordinate, term 0
real(kind_real),intent(in) :: y0      !< Y-coordinate, term 0
real(kind_real),intent(in) :: z0      !< Z-coordinate, term 0
real(kind_real),intent(out) :: output !< Determinant

! Local variables
real(kind_real) :: t1,t2,t3

! Compute determinant terms
t1 = x0*(y1*z2-y2*z1)
t2 = y0*(x1*z2-x2*z1)
t3 = z0*(x1*y2-x2*y1)

! Determinant
output = t1-t2+t3

! Indistinguishability threshold for cross-platform reproducibility
if (small(output,t1).or.small(output,t2).or.small(output,t3)) output = zero

end subroutine stripack_det

!----------------------------------------------------------------------
! Subroutine: stripack_insert
!> Insert K as a neighbor of N1
!----------------------------------------------------------------------
subroutine stripack_insert(k,lp,list,lptr,lnew)

! Passed variables
integer,intent(in) :: k          !< Index of the node to be inserted.
integer,intent(in) :: lp         !< LIST pointer of N2 as a neighbor of N1.
integer,intent(inout) :: list(*) !< Triangulation data structure created by TRMESH. On output,
                                 !< updated with the addition of node K.
integer,intent(inout) :: lptr(*) !< cf. list
integer,intent(inout) :: lnew    !< cf. list

! Local variables
integer :: lsav

lsav = lptr(lp)
lptr(lp) = lnew
list(lnew) = k
lptr(lnew) = lsav
lnew = lnew+1

end subroutine stripack_insert

!----------------------------------------------------------------------
! Subroutine: stripack_intadd
!> Add an interior node to a triangulation
!----------------------------------------------------------------------
subroutine stripack_intadd(kk,i1,i2,i3,list,lptr,lend,lnew)

! Passed variables
integer,intent(in) :: kk         !< Index of the node to be inserted. 1<=KK and KK must not be equal
                                 !< to I1, I2, or I3.
integer,intent(in) :: i1         !< First index of the counterclockwise-ordered sequence of vertices
                                 !< of a triangle which contains node KK.
integer,intent(in) :: i2         !< Second index.
integer,intent(in) :: i3         !< Third index.
integer,intent(inout) :: list(*) !< Triangulation data structure created by TRMESH. Triangle
                                 !< (I1,I2,I3) must be included in the triangulation. On output,
                                 !< updated with the addition of node KK. KK will be connected to
                                 !< nodes I1, I2, and I3.
integer,intent(inout) :: lptr(*) !< cf. list
integer,intent(inout) :: lend(*) !< cf. list
integer,intent(inout) :: lnew    !< cf. list

! Local variables
integer :: k,lp,n1,n2,n3

! Initialization.
k = kk
n1 = i1
n2 = i2
n3 = i3

! Add K as a neighbor of I1, I2, and I3.
call lstptr(lend(n1),n2,list,lptr,lp)
call insert(k,lp,list,lptr,lnew)
call lstptr(lend(n2),n3,list,lptr,lp)
call insert(k,lp,list,lptr,lnew)
call lstptr(lend(n3),n1,list,lptr,lp)
call insert(k,lp,list,lptr,lnew)

! Add I1, I2, and I3 as neighbors of K.
list(lnew) = n1
list(lnew+1) = n2
list(lnew+2) = n3
lptr(lnew) = lnew+1
lptr(lnew+1) = lnew+2
lptr(lnew+2) = lnew
lend(k) = lnew+2
lnew = lnew+3

end subroutine stripack_intadd

!----------------------------------------------------------------------
! Subroutine: stripack_jrand
!> Return a random integer between 1 and N
!----------------------------------------------------------------------
subroutine stripack_jrand(n,ix,iy,iz,output)

! Passed variables
integer,intent(in) :: n       !< Maximum value to be returned.
integer,intent(inout) :: ix   !< First seed initialized to values in the range 1 to 30,000 before
                              !< the first call to JRAND, and not altered between subsequent calls
                              !< (unless a sequence of random numbers is to be repeated by
                              !< reinitializing the seeds).
integer,intent(inout) :: iy   !< Second seed.
integer,intent(inout) :: iz   !< Third seed.
integer,intent(out) :: output !< A random integer in the range 1 to N.

! Local variables
real(kind_real) :: u,x

! Initialization
ix = mod(171*ix,30269)
iy = mod(172*iy,30307)
iz = mod(170*iz,30323)

! Get x
x = (real(ix,kind_real)/30269.0_kind_real) &
 & + (real(iy,kind_real)/30307.0_kind_real) &
 & + (real(iz,kind_real)/30323.0_kind_real)

! Get u
u = x-int(x)

! Random integer
output = int(real(n,kind_real)*u)+1

end subroutine stripack_jrand

!----------------------------------------------------------------------
! Subroutine: stripack_left
!> Determine whether a node is to the left of a plane through the origin
!----------------------------------------------------------------------
subroutine stripack_left(x1,y1,z1,x2,y2,z2,x0,y0,z0,output)

! Passed variables
real(kind_real),intent(in) :: x1 !< X-coordinate, term 1
real(kind_real),intent(in) :: y1 !< Y-coordinate, term 1
real(kind_real),intent(in) :: z1 !< Z-coordinate, term 1
real(kind_real),intent(in) :: x2 !< X-coordinate, term 2
real(kind_real),intent(in) :: y2 !< Y-coordinate, term 2
real(kind_real),intent(in) :: z2 !< Z-coordinate, term 2
real(kind_real),intent(in) :: x0 !< X-coordinate, term 0
real(kind_real),intent(in) :: y0 !< Y-coordinate, term 0
real(kind_real),intent(in) :: z0 !< Z-coordinate, term 0
logical,intent(out) :: output    !< TRUE if and only if N0 is in the closed left hemisphere.

! Local variables
real(kind_real) :: zz

! LEFT = TRUE iff <N0,N1 X N2> = det(N0,N1,N2) >= 0.
call det(x1,y1,z1,x2,y2,z2,x0,y0,z0,zz)
output = sup(zz,zero)

end subroutine stripack_left

!----------------------------------------------------------------------
! Subroutine: stripack_lstptr
!> Return the index of NB in the adjacency list
!----------------------------------------------------------------------
subroutine stripack_lstptr(lpl,nb,list,lptr,output)

! Passed variables
integer,intent(in) :: lpl     !< Equal to LEND(N0).
integer,intent(in) :: nb      !< Index of the node whose pointer is to be returned. NB must be
                              !< connected to N0.
integer,intent(in) :: list(*) !< Triangulation data structure created by TRMESH. Triangle (I1,I2,I3)
                              !< must be included in the triangulation. On output, updated with the
                              !< addition of node KK. KK will be connected to nodes I1, I2, and I3.
integer,intent(in) :: lptr(*) !< cf. list
integer,intent(out) :: output !< Pointer such that LIST(output) = NB or LIST(output) = -NB, unless
                              !< NB is not a neighbor of N0, in which case output = LPL.

! Local variables
integer :: lp,nd

! Initialization
lp = lptr(lpl)

do
   nd = list(lp)
   if (nd==nb) exit
   lp = lptr(lp)
   if (lp==lpl) exit
end do
output = lp

end subroutine stripack_lstptr

!----------------------------------------------------------------------
! Subroutine: stripack_swap
!> Replace the diagonal arc of a quadrilateral with the other diagonal
!----------------------------------------------------------------------
subroutine stripack_swap(in1,in2,io1,io2,list,lptr,lend,lp21)

! Passed variables
integer,intent(in) :: in1        !< First nodal index of the vertices of the quadrilateral. IO1-IO2
                                 !< is replaced by IN1-IN2. (IO1,IO2,IN1) and (IO2,IO1,IN2) must be
                                 !< triangles on input.
integer,intent(in) :: in2        !< Second nodal index.
integer,intent(in) :: io1        !< Third nodal index.
integer,intent(in) :: io2        !< Fourth nodal index.
integer,intent(inout) :: list(*) !< Triangulation data structure created by TRMESH. On output,
                                 !< updated with the swap; triangles (IO1,IO2,IN1) an (IO2,IO1,IN2)
                                 !< are replaced by (IN1,IN2,IO2) and (IN2,IN1,IO1) unless LP21 = 0.
integer,intent(inout) :: lptr(*) !< cf. list
integer,intent(inout) :: lend(*) !< cf. list
integer,intent(out) :: lp21      !< Index of IN1 as a neighbor of IN2 after the swap is performed
                                 !< unless IN1 and IN2 are adjacent on input, in which case
                                 !< LP21 = 0.

! Local variables
integer :: lp,lph,lpsav

! Test for IN1 and IN2 adjacent.
call lstptr(lend(in1),in2,list,lptr,lp)
if (abs(list(lp))==in2) then
   lp21 = 0
   return
end if

! Delete IO2 as a neighbor of IO1.
call lstptr(lend(io1),in2,list,lptr,lp)
lph = lptr(lp)
lptr(lp) = lptr(lph)

! If IO2 is the last neighbor of IO1, make IN2 the last neighbor.
if (lend(io1)==lph) lend(io1) = lp

! Insert IN2 as a neighbor of IN1 following IO1 using the hole created above.
call lstptr(lend(in1),io1,list,lptr,lp)
lpsav = lptr(lp)
lptr(lp) = lph
list(lph) = in2
lptr(lph) = lpsav

! Delete IO1 as a neighbor of IO2.
call lstptr(lend(io2),in1,list,lptr,lp)
lph = lptr(lp)
lptr(lp) = lptr(lph)

! If IO1 is the last neighbor of IO2, make IN1 the last neighbor.
if (lend(io2)==lph) lend(io2) = lp

! Insert IN1 as a neighbor of IN2 following IO2.
call lstptr(lend(in2),io2,list,lptr,lp)
lpsav = lptr(lp)
lptr(lp) = lph
list(lph) = in1
lptr(lph) = lpsav
lp21 = lph

end subroutine stripack_swap

!----------------------------------------------------------------------
! Subroutine: stripack_swptst
!> Decide whether to replace a diagonal arc by the other
!----------------------------------------------------------------------
subroutine stripack_swptst(n1,n2,n3,n4,x,y,z,output)

integer,intent(in) :: n1           !< First index of the four nodes defining the quadrilateral with
                                   !< N1 adjacent to N2, and (N1,N2,N3) in counterclockwise order.
                                   !< The arc connecting N1 to N2 should be replaced by an arc
                                   !< connecting N3 to N4 if SWPTST = TRUE. Refer to subroutine
                                   !< SWAP.
integer,intent(in) :: n2           !< Second index.
integer,intent(in) :: n3           !< Third index.
integer,intent(in) :: n4           !< Fourth index.
real(kind_real),intent(in) :: x(*) !< X-coordinate of the nodes.
real(kind_real),intent(in) :: y(*) !< Y-coordinate of the nodes.
real(kind_real),intent(in) :: z(*) !< Z-coordinate of the nodes.
logical,intent(out) :: output      !< TRUE if and only if the arc connecting N1 and N2 should be
                                   !< swapped for an arc connecting N3 and N4.

! Local variables
real(kind_real) :: dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,x4,y4,z4,zz

! Initialization
x4 = x(n4)
y4 = y(n4)
z4 = z(n4)
dx1 = x(n1)-x4
dx2 = x(n2)-x4
dx3 = x(n3)-x4
dy1 = y(n1)-y4
dy2 = y(n2)-y4
dy3 = y(n3)-y4
dz1 = z(n1)-z4
dz2 = z(n2)-z4
dz3 = z(n3)-z4

! N4 lies above the plane of (N1,N2,N3) iff N3 lies above the plane of (N2,N1,N4) iff
! Det(N3-N4,N2-N4,N1-N4) = (N3-N4,N2-N4 X N1-N4) > 0.
call det(dx2,dy2,dz2,dx1,dy1,dz1,dx3,dy3,dz3,zz)
output = sup(zz,zero)

end subroutine stripack_swptst

!----------------------------------------------------------------------
! Subroutine: stripack_trfind
!> Locate a point relative to a triangulation
!----------------------------------------------------------------------
subroutine stripack_trfind(nst,p,n,x,y,z,list,lptr,lend,b1,b2,b3,i1,i2,i3)

integer,intent(in) :: nst          !< Index of a node at which TRFIND begins its search. Search time
                                   !< depends on the proximity of this node to P.
real(kind_real) :: p(3)            !< X, Y, and Z coordinates (in that order) of the point P to be
                                   !< located.
integer,intent(in) :: n            !< Number of nodes in the triangulation 3<=N.
real(kind_real),intent(in) :: x(n) !< X-coordinate of the triangulation nodes (unit vector).
real(kind_real),intent(in) :: y(n) !< Y-coordinate of the triangulation nodes (unit vector).
real(kind_real),intent(in) :: z(n) !< Z-coordinate of the triangulation nodes (unit vector).
integer,intent(in) :: list(*)      !< Triangulation data structure created by TRMESH.
integer,intent(in) :: lptr(*)      !< cf. list
integer,intent(in) :: lend(*)      !< cf. list
real(kind_real),intent(out) :: b1  !< First unnormalized barycentric coordinate of the central
                                   !< projection of P onto the underlying planar triangle if P is in
                                   !< the convex hull of the nodes. These parameters are not altered
                                   !< if I1 = 0.
real(kind_real),intent(out) :: b2  !< Second unnormalized barycentric coordinate.
real(kind_real),intent(out) :: b3  !< Third unnormalized barycentric coordinate.
integer,intent(out) :: i1          !< First counterclockwise-ordered vertex index of a triangle
                                   !< containing P if P is contained in a triangle. If P is not in
                                   !< the convex hull of the nodes, I1 and I2 are the rightmost and
                                   !< leftmost (boundary) nodes that are visible from P, and I3 = 0.
                                   !< (If all boundary nodes are visible from P, then I1 and I2
                                   !< coincide.) I1 = I2 = I3 = 0 if P and all of the nodes are
                                   !< coplanar (lie on a common great circle).
integer,intent(out) :: i2          !< Second counterclockwise-ordered vertex index.
integer,intent(out) :: i3          !< Third counterclockwise-ordered vertex index.

! Local variables
integer :: lp,n0,n1,n1s,n2,n2s,n3,n4,next,nf,nl
integer,save :: ix = 1,iy = 2,iz = 3
real(kind_real) :: output,eps,ptn1,ptn2,q(3),s12,tol,xp,yp,zp

! Initialize variables.
xp = p(1)
yp = p(2)
zp = p(3)
n0 = nst
if ((n0<1).or.(n<n0)) call jrand(n,ix,iy,iz,n0)

! Compute the relative machine precision EPS and TOL.
eps = epsilon(eps)
tol = hundred*eps

! Set NF and NL to the first and last neighbors of N0, and initialize N1 = NF.
2 continue
lp = lend(n0)
nl = list(lp)
lp = lptr(lp)
nf = list(lp)
n1 = nf

! Find a pair of adjacent neighbors N1,N2 of N0 that define a wedge containing P:  P LEFT N0->N1 and
! P RIGHT N0->N2.
if (0<nl) then
   ! N0 is an interior node. Find N1.
   3 continue
   call det(x(n0),y(n0),z(n0),x(n1),y(n1),z(n1),xp,yp,zp,output)
   if (inf(output,zero)) then
      lp = lptr(lp)
      n1 = list(lp)
      if (n1 == nl) go to 6
      go to 3
   end if
else
   ! N0 is a boundary node. Test for P exterior.
   nl = -nl

   ! Is P to the right of the boundary edge N0->NF?
   call det(x(n0),y(n0),z(n0),x(nf),y(nf),z(nf),xp,yp,zp,output)
   if (inf(output,zero)) then
      n1 = n0
      n2 = nf
      go to 9
   end if

   ! Is P to the right of the boundary edge NL->N0?
   call det(x(nl),y(nl),z(nl),x(n0),y(n0),z(n0),xp,yp,zp,output)
   if (inf(output,zero)) then
      n1 = nl
      n2 = n0
      go to 9
   end if
end if

! P is to the left of arcs N0->N1 and NL->N0. Set N2 to the next neighbor of N0 (following N1).
4 continue
lp = lptr(lp)
n2 = abs(list(lp))
call det(x(n0),y(n0),z(n0),x(n2),y(n2),z(n2),xp,yp,zp,output)
if (inf(output,zero)) go to 7
n1 = n2
if (n1/=nl) go to 4
call det(x(n0),y(n0),z(n0),x(nf),y(nf),z(nf),xp,yp,zp,output)
if (inf(output,zero)) go to 6

! P is left of or on arcs N0->NB for all neighbors NB of N0. Test for P = +/-N0.
if (inf(abs(x(n0)*xp+y(n0)*yp+z(n0)*zp),one-four*eps)) then
   ! All points are collinear iff P Left NB->N0 for all neighbors NB of N0. Search the neighbors of
   ! N0. Note:  N1 = NL and LP points to NL.
   do
      call det(x(n1),y(n1),z(n1),x(n0),y(n0),z(n0),xp,yp,zp,output)
      if (inf(output,zero)) exit
      lp = lptr(lp)
      n1 = abs(list(lp))
      if (n1==nl) then
         i1 = 0
         i2 = 0
         i3 = 0
         return
      end if
   end do
end if

! P is to the right of N1->N0, or P = +/-N0. Set N0 to N1 and start over.
n0 = n1
go to 2

! P is between arcs N0->N1 and N0->NF.
6 continue
n2 = nf

! P is contained in a wedge defined by geodesics N0-N1 and N0-N2, where N1 is adjacent to N2. Save
! N1 and N2 to test for cycling.
7 continue
n3 = n0
n1s = n1
n2s = n2

! Top of edge-hopping loop:
8 continue
call det(x(n1),y(n1),z(n1),x(n2),y(n2),z(n2),xp,yp,zp,b3)
if (inf(b3,zero)) then
   ! Set N4 to the first neighbor of N2 following N1 (the node opposite N2->N1) unless N1->N2 is a
   ! boundary arc.
   call lstptr(lend(n2),n1,list,lptr,lp)
   if (list(lp)<0) go to 9
   lp = lptr(lp)
   n4 = abs(list(lp))

   !  Define a new arc N1->N2 which intersects the geodesic N0-P.
   call det(x(n0),y(n0),z(n0),x(n4),y(n4),z(n4),xp,yp,zp,output)
   if (inf(output,zero)) then
      n3 = n2
      n2 = n4
      n1s = n1
      if ((n2/=n2s).and.(n2/=n0)) go to 8
   else
      n3 = n1
      n1 = n4
      n2s = n2
      if ((n1/=n1s).and.(n1/=n0)) go to 8
   end if

   ! The starting node N0 or edge N1-N2 was encountered again, implying a cycle (infinite loop).
   ! Restart with N0 randomly selected.
   call jrand(n,ix,iy,iz,n0)
   go to 2
end if

! P is in (N1,N2,N3) unless N0, N1, N2, and P are collinear or P is close to -N0.
if (supeq(b3,eps)) then
   ! B3 /= 0.
    call det(x(n2),y(n2),z(n2),x(n3),y(n3),z(n3),xp,yp,zp,b1)
    call det(x(n3),y(n3),z(n3),x(n1),y(n1),z(n1),xp,yp,zp,b2)

   ! Restart with N0 randomly selected.
   if (inf(b1,-tol).or.inf(b2,-tol)) then
      call jrand(n,ix,iy,iz,n0)
      go to 2
   end if
else
    ! B3 = 0 and thus P lies on N1->N2. Compute B1 = Det(P,N2 X N1,N2) and B2 = Det(P,N1,N2 X N1).
    b3 = zero
    s12 = x(n1)*x(n2)+y(n1)*y(n2)+z(n1)*z(n2)
    ptn1 = xp*x(n1)+yp*y(n1)+zp*z(n1)
    ptn2 = xp*x(n2)+yp*y(n2)+zp*z(n2)
    b1 = ptn1-s12*ptn2
    b2 = ptn2-s12*ptn1

   ! Restart with N0 randomly selected.
   if (inf(b1,-tol).or.inf(b2,-tol)) then
      call jrand(n,ix,iy,iz,n0)
      go to 2
   end if
end if

! P is in (N1,N2,N3).
i1 = n1
i2 = n2
i3 = n3
b1 = max(b1,zero)
b2 = max(b2,zero)
return

! P Right N1->N2, where N1->N2 is a boundary edge. Save N1 and N2, and set NL = 0 to indicate that
! NL has not yet been found.
9 continue
n1s = n1
n2s = n2
nl = 0

! Counterclockwise Boundary Traversal:
10 continue
lp = lend(n2)
lp = lptr(lp)
next = list(lp)
call det(x(n2),y(n2),z(n2),x(next),y(next),z(next),xp,yp,zp,output)
if (supeq(output,zero)) then
   ! N2 is the rightmost visible node if P Forward N2->N1 or NEXT Forward N2->N1.
   ! Set Q to (N2 X N1) X N2.
    s12 = x(n1)*x(n2)+y(n1)*y(n2)+z(n1)*z(n2)
    q(1) = x(n1)-s12*x(n2)
    q(2) = y(n1)-s12*y(n2)
    q(3) = z(n1)-s12*z(n2)
    if (sup(xp*q(1)+yp*q(2)+zp*q(3),zero).or.sup(x(next)*q(1)+y(next)*q(2)+z(next)*q(3),zero)) &
      go to 11

   ! N1, N2, NEXT, and P are nearly collinear, and N2 is the leftmost visible node.
   nl = n2
end if

! Bottom of counterclockwise loop:
n1 = n2
n2 = next
if (n2/=n1s) go to 10

! All boundary nodes are visible from P.
i1 = n1s
i2 = n1s
i3 = 0
return

! N2 is the rightmost visible node.
11 continue
nf = n2
if (nl==0) then
   ! Restore initial values of N1 and N2, and begin the search for the leftmost visible node.
   n2 = n2s
   n1 = n1s

   ! Clockwise Boundary Traversal:
   12  continue
   lp = lend(n1)
   next = -list(lp)
   call det(x(next),y(next),z(next),x(n1),y(n1),z(n1),xp,yp,zp,output)
   if (supeq(output,zero)) then
      ! N1 is the leftmost visible node if P or NEXT is forward of N1->N2.
      ! Compute Q = N1 X (N2 X N1).
      s12 = x(n1)*x(n2)+y(n1)*y(n2)+z(n1)*z(n2)
      q(1) = x(n2)-s12*x(n1)
      q(2) = y(n2)-s12*y(n1)
      q(3) = z(n2)-s12*z(n1)
      if (sup(xp*q(1)+yp*q(2)+zp*q(3),zero).or.sup(x(next)*q(1)+y(next)*q(2)+z(next)*q(3),zero)) &
        go to 13

      ! P, NEXT, N1, and N2 are nearly collinear and N1 is the rightmost visible node.
      nf = n1
   end if

   ! Bottom of clockwise loop:
   n2 = n1
   n1 = next
   if (n1/=n1s) go to 12

   ! All boundary nodes are visible from P.
   i1 = n1
   i2 = n1
   i3 = 0
   return

   ! N1 is the leftmost visible node.
   13  continue
   nl = n1
end if

! NF and NL have been found.
i1 = nf
i2 = nl
i3 = 0

end subroutine stripack_trfind

!----------------------------------------------------------------------
! Subroutine: stripack_trmesh
!> Create a Delaunay triangulation on the unit sphere
!----------------------------------------------------------------------
subroutine stripack_trmesh(n,x,y,z,list,lptr,lend,lnew,near,next,dist,ier)

integer,intent(in) :: n                  !< Number of nodes in the triangulation. 3<=N.
real(kind_real),intent(in) :: x(n)       !< X-coordinate of distinct nodes. (X(K),Y(K), Z(K)) is
                                         !< referred to as node K, and K is referred to as a nodal
                                         !< index. It is required that X(K)**2+Y(K)**2+Z(K)**2 = 1
                                         !< for all K. The first three nodes must not be colinear
                                         !< (lie on a common great circle).
real(kind_real),intent(in) :: y(n)       !< Y-coordinate of distinct nodes.
real(kind_real),intent(in) :: z(n)       !< Z-coordinate of distinct nodes.
integer,intent(out) :: list(6*(n-2))     !< Nodal indexes which, along with LPTR, LEND, and LNEW,
                                         !< define the triangulation as a set of N adjacency lists;
                                         !< counterclockwise-ordered sequences of neighboring nodes
                                         !< such that the first and last neighbors of a boundary
                                         !< node are boundary nodes (the first neighbor of an
                                         !< interior node is arbitrary). In order to distinguish
                                         !< between interior and boundary nodes, the last neighbor
                                         !< of each boundary node is represented by the negative of
                                         !< its index.
integer,intent(out) :: lptr(6*(n-2))     !< Set of pointers (LIST indexes) in one-to-one
                                         !< correspondence with the elements of LIST. LIST(LPTR(I))
                                         !< indexes the node which follows LIST(I) in cyclical
                                         !< counterclockwise order (the first neighbor follows the
                                         !< last neighbor).
integer,intent(out) :: lend(n)           !< Pointers to adjacency lists. LEND(K) points to the last
                                         !< neighbor of node K. LIST(LEND(K))<0 if and only if K is
                                         !< a boundary node.
integer,intent(out) :: lnew              !< Pointer to the first empty location in LIST and LPTR
                                         !< (list length plus one).
integer,intent(inout) :: near(n)         !< Workspace used to efficiently determine the nearest
                                         !< triangulation node to each unprocessed node for use by
                                         !< ADDNOD.
integer,intent(inout) :: next(n)         !< Workspace used to efficiently determine the nearest
                                         !< triangulation node to each unprocessed node for use by
                                         !< ADDNOD.
real(kind_real),intent(inout) :: dist(n) !< Workspace used to efficiently determine the nearest
                                         !< triangulation node to each unprocessed node for use by
                                         !< ADDNOD.
integer,intent(out) :: ier               !< error return code
                                         !<     0, if no errors were encountered.
                                         !<    -1, if N < 3 on input.
                                         !<    -2, if the first three nodes are collinear.
                                         !<     L, if nodes L and M coincide for some L < M.  The data structure 
                                         !<      represents a triangulation of nodes 1 to M-1 in this case.


! Local variables
integer :: i,i0,j,k,lp,lpl,nexti,nn
real(kind_real) :: d,d1,d2,d3
logical :: output_1,output_2

! Initialization
ier = 0
nn = n
if (nn<3) then
   call abor1_ftn('stripack_trmesh: N < 3')
   ier = -1;
   return
endif
list = 0
lptr = 0
lend = 0

! Store the first triangle in the linked list.
call left(x(1),y(1),z(1),x(2),y(2),z(2),x(3),y(3),z(3),output_1)
call left(x(2),y(2),z(2),x(1),y(1),z(1),x(3),y(3),z(3),output_2)

if (.not.output_1) then
   ! The first triangle is (3,2,1) = (2,1,3) = (1,3,2).
   list(1) = 3
   lptr(1) = 2
   list(2) = -2
   lptr(2) = 1
   lend(1) = 2

   list(3) = 1
   lptr(3) = 4
   list(4) = -3
   lptr(4) = 3
   lend(2) = 4

   list(5) = 2
   lptr(5) = 6
   list(6) = -1
   lptr(6) = 5
   lend(3) = 6
else if (.not.output_2) then
   ! The first triangle is (1,2,3):  3 Strictly Left 1->2, i.e., node 3 lies in the left hemisphere
   ! defined by arc 1->2.
   list(1) = 2
   lptr(1) = 2
   list(2) = -3
   lptr(2) = 1
   lend(1) = 2

   list(3) = 3
   lptr(3) = 4
   list(4) = -1
   lptr(4) = 3
   lend(2) = 4

   list(5) = 1
   lptr(5) = 6
   list(6) = -2
   lptr(6) = 5
   lend(3) = 6
else
   ! The first three nodes are colinear.
   call abor1_ftn('stripack_trmesh: The first three nodes are colinear')
   ier = -2
   return
end if

! Initialize LNEW and test for N = 3.
lnew = 7
if (nn==3) then
   ier = 0
   return
end if

! A nearest-node data structure (NEAR, NEXT, and DIST) is used to obtain an expected-time
! (N*log(N)) incremental algorithm by enabling constant search time for locating each new node in
! the triangulation. For each unprocessed node K, NEAR(K) is the index of the triangulation node
! closest to K (used as the starting point for the search in Subroutine TRFIND) and DIST(K) is an
! increasing function of the arc length (angular distance) between nodes K and NEAR(K): -cos(a) for
! arc length a. Since it is necessary to efficiently find the subset of unprocessed nodes associated
! with each triangulation node J (those that have J as their NEAR entries), the subsets are stored
! in NEAR and NEXT as follows:  for each node J in the triangulation, I = NEAR(J) is the first
! unprocessed node in J's set (with I = 0 if the set is empty), L = NEXT(I) (if 0<I) is the second,
! NEXT(L) (if 0<L) is the third, etc. The nodes in each set are initially ordered by increasing
! indexes (which maximizes efficiency) but that ordering is not maintained as the data structure is
! updated.

! Initialize the data structure for the single triangle.
near(1) = 0
near(2) = 0
near(3) = 0
do k=nn,4,-1
   d1 = -(x(k)*x(1)+y(k)*y(1)+z(k)*z(1))
   d2 = -(x(k)*x(2)+y(k)*y(2)+z(k)*z(2))
   d3 = -(x(k)*x(3)+y(k)*y(3)+z(k)*z(3))
   if (inf(d1,d2).and.inf(d1,d3)) then
      near(k) = 1
      dist(k) = d1
      next(k) = near(1)
      near(1) = k
   else if (inf(d2,d1).and.inf(d2,d3)) then
      near(k) = 2
      dist(k) = d2
      next(k) = near(2)
      near(2) = k
   else
      near(k) = 3
      dist(k) = d3
      next(k) = near(3)
      near(3) = k
   end if
end do

! Add the remaining nodes.
do k=4,nn
   call addnod(near(k),k,x,y,z,list,lptr,lend,lnew,ier)
   if (ier/=0) then
      return
   endif

   ! Remove K from the set of unprocessed nodes associated with NEAR(K).
   i = near(k)
   i0 = i
   if (near(i)==k) then
      near(i) = next(k)
   else
      i = near(i)
      do
         i0 = i
         i = next(i0)
         if (i==k) exit
      end do
      next(i0) = next(k)
   end if
   near(k) = 0

   ! Loop on neighbors J of node K.
   lpl = lend(k)
   lp = lpl
   3 continue
   lp = lptr(lp)
   j = abs(list(lp))

   ! Loop on elements I in the sequence of unprocessed nodes associated with J:  K is a candidate
   ! for replacing J as the nearest triangulation node to I. The next value of I in the sequence,
   ! NEXT(I), must be saved before I is moved because it is altered by adding I to K's set.
   i = near(j)
   do
      if (i==0) exit
      nexti = next(i)

      ! Test for the distance from I to K less than the distance from I to J. Indistinguishability
      ! threshold for cross-platform reproducibility
      d =-(x(i)*x(k)+y(i)*y(k)+z(i)*z(k))
      if (inf(d,dist(i))) then
         ! Replace J by K as the nearest triangulation node to I: update NEAR(I) and DIST(I), and
         ! remove I from J's set of unprocessed nodes and add it to K's set.
         near(i) = k
         dist(i) = d
         if (i==near(j)) then
            near(j) = nexti
         else
            next(i0) = nexti
         end if
         next(i) = near(k)
         near(k) = i
      else
         i0 = i
      end if
      i = nexti
   end do

   ! Bottom of loop on neighbors J.
   if (lp/=lpl) go to 3
end do

end subroutine stripack_trmesh




subroutine stripack_trlist2 ( n, list, lptr, lend, nt, ltri, ier )

   !*****************************************************************************80
   !
   !! TRLIST2 converts a triangulation data structure to a triangle list.
   !
   !  Discussion:
   !
   !    This subroutine converts a triangulation data structure
   !    from the linked list created by TRMESH to a triangle list.
   !
   !    It is a version of TRLIST for the special case where the triangle
   !    list should only include the nodes that define each triangle.
   !
   !  Modified:
   !
   !    21 July 2007
   !
   !  Author:
   !
   !    Robert Renka
   !
   !  Reference:
   !
   !    Robert Renka,
   !    Algorithm 772: STRIPACK,
   !    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
   !    ACM Transactions on Mathematical Software,
   !    Volume 23, Number 3, September 1997, pages 416-434.
   !
   !  Parameters:
   !
   !    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
   !    3 <= N.
   !
   !    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), linked
   !    list data structure defining the triangulation.  Refer to TRMESH.
   !
   !    Output, integer ( kind = 4 ) NT, the number of triangles in the
   !    triangulation unless IER /=0, in which case NT = 0.  NT = 2N-NB-2 if
   !    NB >= 3 or 2N-4 if NB = 0, where NB is the number of boundary nodes.
   !
   !    Output, integer ( kind = 4 ) LTRI(3,*).  The second dimension of LTRI
   !    must be at least NT, where NT will be at most 2*N-4.  The J-th column
   !    contains the vertex nodal indexes associated with triangle J for
   !    J = 1,...,NT.  The vertices are ordered counterclockwise with the first
   !    vertex taken to be the one with smallest index.  Thus, LTRI(2,J) and
   !    LTRI(3,J) are larger than LTRI(1,J) and index adjacent neighbors of node
   !    LTRI(1,J).  The triangles are ordered on first (smallest) vertex indexes.
   !
   !    Output, integer ( kind = 4 ) IER, error indicator.
   !    0, if no errors were encountered.
   !    1, if N is outside its valid range on input.
   !    2, if the triangulation data structure (LIST,LPTR,LEND) is invalid.
   !      Note, however, that these arrays are not completely tested for validity.
   !
   !  Local parameters:
   !
   !    I,J =      LTRI row indexes (1 to 3) associated with
   !               triangles KT and KN, respectively
   !    I1,I2,I3 = Nodal indexes of triangle KN
   !    ISV =      Variable used to permute indexes I1,I2,I3
   !    KA =       Arc index and number of currently stored arcs
   !    KN =       Index of the triangle that shares arc I1-I2 with KT
   !    KT =       Triangle index and number of currently stored triangles
   !    LP =       LIST pointer
   !    LP2 =      Pointer to N2 as a neighbor of N1
   !    LPL =      Pointer to the last neighbor of I1
   !    LPLN1 =    Pointer to the last neighbor of N1
   !    N1,N2,N3 = Nodal indexes of triangle KT
   !    NM2 =      N-2
   !
       implicit none
   
       integer ( kind = 4 ) n
       integer ( kind = 4 ) i
       integer ( kind = 4 ) i1
       integer ( kind = 4 ) i2
       integer ( kind = 4 ) i3
       integer ( kind = 4 ) ier
       integer ( kind = 4 ) isv
       integer ( kind = 4 ) j
       integer ( kind = 4 ) ka
       integer ( kind = 4 ) kn
       integer ( kind = 4 ) kt
       integer ( kind = 4 ) lend(n)
       integer ( kind = 4 ) list(6*(n-2))
       integer ( kind = 4 ) lp
       integer ( kind = 4 ) lp2
       integer ( kind = 4 ) lpl
       integer ( kind = 4 ) lpln1
       integer ( kind = 4 ) lptr(6*(n-2))
       integer ( kind = 4 ) ltri(3,2*n-4)
       integer ( kind = 4 ) n1
       integer ( kind = 4 ) n2
       integer ( kind = 4 ) n3
       integer ( kind = 4 ) nm2
       integer ( kind = 4 ) nt
   !
   !  Test for invalid input parameters.
   !
       if ( n < 3 ) then
         nt = 0
         ier = 1
         return
       end if
   !
   !  Initialize parameters for loop on triangles KT = (N1,N2,
   !  N3), where N1 < N2 and N1 < N3.
   !
   !  KA,KT = Numbers of currently stored arcs and triangles.
   !  NM2 = Upper bound on candidates for N1.
   !
       ka = 0
       kt = 0
       nm2 = n-2
   !
   !  Loop on nodes N1.
   !
       do n1 = 1, nm2
   !
   !  Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
   !  to the last neighbor of N1, and LP2 points to N2.
   !
         lpln1 = lend(n1)
         lp2 = lpln1
   
   1   continue
   
       lp2 = lptr(lp2)
       n2 = list(lp2)
       lp = lptr(lp2)
       n3 = abs ( list(lp) )
   
       if ( n2 < n1 .or. n3 < n1 ) then
         go to 8
       end if
   
   !
   !  Add a new triangle KT = (N1,N2,N3).
   !
       kt = kt + 1
       ltri(1,kt) = n1
       ltri(2,kt) = n2
       ltri(3,kt) = n3
   !
   !  Loop on triangle sides (I2,I1) with neighboring triangles
   !  KN = (I1,I2,I3).
   !
       do i = 1, 3
   
         if ( i == 1 ) then
           i1 = n3
           i2 = n2
         else if ( i == 2 ) then
           i1 = n1
           i2 = n3
         else
           i1 = n2
           i2 = n1
         end if
   !
   !  Set I3 to the neighbor of I1 that follows I2 unless
   !  I2->I1 is a boundary arc.
   !
         lpl = lend(i1)
         lp = lptr(lpl)
   
         do
   
           if ( list(lp) == i2 ) then
             go to 3
           end if
   
           lp = lptr(lp)
   
           if ( lp == lpl ) then
             exit
           end if
   
         end do
   !
   !  Invalid triangulation data structure:  I1 is a neighbor of
   !  I2, but I2 is not a neighbor of I1.
   !
         if ( abs ( list(lp) ) /= i2 ) then
           nt = 0
           ier = 2
           return
         end if
   !
   !  I2 is the last neighbor of I1.  Bypass the search for a neighboring
   !  triangle if I2->I1 is a boundary arc.
   !
         kn = 0
   
         if ( list(lp) < 0 ) then
           go to 6
         end if
   !
   !  I2->I1 is not a boundary arc, and LP points to I2 as
   !  a neighbor of I1.
   !
   3       continue
   
         lp = lptr(lp)
         i3 = abs ( list(lp) )
   !
   !  Find J such that LTRI(J,KN) = I3 (not used if KT < KN),
   !  and permute the vertex indexes of KN so that I1 is smallest.
   !
         if ( i1 < i2 .and. i1 < i3 ) then
           j = 3
         else if ( i2 < i3 ) then
           j = 2
           isv = i1
           i1 = i2
           i2 = i3
           i3 = isv
         else
           j = 1
           isv = i1
           i1 = i3
           i3 = i2
           i2 = isv
         end if
   !
   !  Test for KT < KN (triangle index not yet assigned).
   !
         if ( n1 < i1 ) then
           cycle
         end if
   !
   !  Find KN, if it exists, by searching the triangle list in
   !  reverse order.
   !
         do kn = kt-1, 1, -1
           if ( ltri(1,kn) == i1 .and. &
                ltri(2,kn) == i2 .and. &
                ltri(3,kn) == i3 ) then
             go to 5
           end if
         end do
   
         cycle
   
   5       continue
   
   6       continue
   
         end do
   !
   !  Bottom of loop on triangles.
   !
   8     continue
   
       if ( lp2 /= lpln1 ) then
         go to 1
       end if
   
   9     continue
   
       end do
   
       nt = kt
       ier = 0
   
       return
       end



end module stripack_mod
