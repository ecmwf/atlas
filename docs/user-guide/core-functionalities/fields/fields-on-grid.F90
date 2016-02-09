program main

use atlas_module

integer, parameter      :: wp = 8
integer                 :: jnode
character(len=1024)     :: gridID
type(atlas_ReducedGrid) :: reducedGrid
type(atlas_Field)       :: windfield
real(wp), pointer       :: wind(:,:)

call atlas_init()

call atlas_resource("--grid", "N32", gridID)
reducedGrid = atlas_ReducedGrid(gridID)

windfield = atlas_Field("wind", &
                      & atlas_real(wp), &
                      & [2, 100])

call windfield%data(wind)
do jnode=1,100
  wind(1, jnode) = 0.01_wp + real(jnode,kind=wp)
  wind(2, jnode) = 0.02_wp + real(jnode,kind=wp)
enddo

write(unit=6, fmt = "(A,F16.8,1x)") "wind(1, 10) = ", wind(1,10)
write(unit=6, fmt = "(A,F16.8,1x)") "wind(2, 10) = ", wind(2,10)

call reducedGrid%final()
call windfield%final()

call atlas_finalize()

end program main

