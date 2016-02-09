program main

use atlas_module

! IC parameters
integer, parameter              :: wp = 8
! pi and commonly used divisions
real(wp), parameter             :: rpi = 2.0_wp * asin(1.0_wp)
real(wp), parameter             :: deg2rad = rpi / 180.
real(wp)                        :: zlatc = 0._wp * rpi
real(wp)                        :: zlonc = 1._wp * rpi
real(wp)                        :: zrad  = 2._wp * rpi / 9._wp
real(wp)                        :: zdist, zlon, zlat

integer                         :: jnode
character(len=1024)             :: gridID
type(atlas_ReducedGrid)         :: reducedGrid
type(atlas_Field)               :: pressureField
real(wp), pointer               :: pressure(:)

call atlas_init()

call atlas_resource("--grid", "N32", gridID)
reducedGrid = atlas_ReducedGrid(gridID)

pressureField = atlas_Field("pressure", atlas_real(wp), [reducedGrid%npts()])
call pressureField%data(pressure)

jnode = 1
do jlat=1,reducedGrid%nlat()
  zlat = reducedGrid%lat(jlat)
  zlat = zlat * deg2rad
  do jlon=1,reducedGrid%nlon(jlat)
    zlon = reducedGrid%lon(jlat,jlon)
    zlon = zlon * deg2rad

    zdist = 2._wp * sqrt((cos(zlat) * sin((zlon - zlonc) / 2))**2 + &
                      & sin((zlat - zlatc) / 2)**2)

    pressure(jnode) = 0._wp
    if (zdist < zrad) then
      pressure(jnode) = 0.5_wp * (1._wp + cos(rpi * zdist / zrad))
    endif
    jnode = jnode + 1
  enddo
enddo

write(6, *) "================================================================="
write(6, *) "memory pressureField = ", pressureField%bytes()/1000000000., "GB"
write(6, *) "================================================================="
do jnode=12,18
  write(6, *) "node = ", jnode, "pressure = ", pressure(jnode)
enddo

call reducedGrid  %final()
call pressureField%final()

call atlas_finalize()

end program main

