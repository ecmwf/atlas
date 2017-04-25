program main

use, intrinsic :: iso_c_binding, only : c_double
use atlas_module
implicit none

integer, parameter            :: wp = c_double
integer                       :: jnode
character(len=1024)           :: string
character(len=:), allocatable :: unitsW, unitsP
type(atlas_Field)             :: field_pressure , field_wind
type(atlas_Field)             :: field_pressure2, field_wind2
real(wp), pointer             :: pressure(:), wind    (:,:)
type(atlas_FieldSet)          :: fields
type(atlas_Metadata)          :: metadata
call atlas_library%initialise()

! Define fields
field_pressure = atlas_Field(name="pressure", kind=atlas_real(wp), shape=[100])
field_wind     = atlas_Field(name="wind"    , kind=atlas_real(wp), shape=[2, 100])

! Access fields data
call field_pressure%data(pressure)
call field_wind    %data(wind)

! Assign values to fields
do jnode=1,100
  pressure(jnode) = 101325._wp
  wind(1, jnode)  = 0.01_wp + real(jnode,kind=wp)
  wind(2, jnode)  = 0.02_wp + real(jnode,kind=wp)
enddo

! Add info to fields
metadata = field_pressure%metadata()
call metadata%set("units", "[Pa]")
call metadata%get("units", unitsP)
metadata = field_wind%metadata()
call metadata%set("units", "[m/s]")
call metadata%get("units", unitsW)

! Define fieldSet
fields = atlas_FieldSet()
call fields%add(field_pressure) ! Add field_pressure to fieldSet
call fields%add(field_wind)     ! Add field_wind to fieldSet

! Retrieve field from fieldSet
field_pressure2 = fields%field("pressure")
field_wind2     = fields%field("wind")

! Print some useful info
write(string, *) "name   = ", field_wind%name()
call atlas_log%info(string)
write(string, *) "size   = ", field_wind%size()
call atlas_log%info(string)
write(string, *) "units  = ", unitsW
call atlas_log%info(string)
write(string, *) "rank   = ", field_wind%rank()
call atlas_log%info(string)
write(string, *) "shape(1)  = ", field_wind%shape(1)
call atlas_log%info(string)
write(string, *) "shape(2)  = ", field_wind%shape(2)
call atlas_log%info(string)
write(string, *) "shape  = ", field_wind%shape()
call atlas_log%info(string)
write(string, *) "memory = ", field_wind%bytes(), "bytes"
call atlas_log%info(string)
write(string, *) "type   = ", field_wind%datatype()
call atlas_log%info(string)
write(string, *) "kind   = ", field_wind%kind()
call atlas_log%info(string)

! Print some values
write(string, *) "pressure(10) = ", pressure(10)
call atlas_log%info(string)
write(string, *) "wind(1, 10)  = ", wind(1,10)
call atlas_log%info(string)
write(string, *) "wind(2, 10)  = ", wind(2,10)
call atlas_log%info(string)

! Finalize object to release memory
call field_pressure%final()
call field_wind    %final()
call fields        %final()

call atlas_library%finalise()
end program main
