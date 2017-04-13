program main
use, intrinsic :: iso_c_binding, only : c_double
use atlas_module
implicit none

integer, parameter              :: wp = c_double
real(wp), parameter             :: rpi = 2.0_wp * asin(1.0_wp)
real(wp), parameter             :: deg2rad = rpi / 180.
real(wp)                        :: zlatc = 0._wp * rpi
real(wp)                        :: zlonc = 1._wp * rpi
real(wp)                        :: zrad  = 2._wp * rpi / 9._wp
real(wp)                        :: zdist, zlon, zlat

integer                         :: jnode, jlon, jlat
character(len=1024)             :: string
type(atlas_StructuredGrid)     :: grid
type(atlas_Field)               :: field_pressure
real(wp), pointer               :: pressure(:)

call atlas_library%initialise()

grid = atlas_StructuredGrid( "O32" )

field_pressure = atlas_Field("pressure", atlas_real(wp), [grid%size()])
call field_pressure%data(pressure)

jnode = 1
do jlat=1,grid%ny()
  zlat = grid%y(jlat)
  zlat = zlat * deg2rad
  do jlon=1,grid%nx(jlat)
    zlon = grid%x(jlon,jlat)
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

write(string, *) "===================================================="
call atlas_log%info(string)
write(string, *) "memory field_pressure = ", &
            & field_pressure%bytes()/1000000000., "GB"
call atlas_log%info(string)
write(string, *) "===================================================="
call atlas_log%info(string)

call grid%final()
call field_pressure%final()

call atlas_library%finalise()

end program main

