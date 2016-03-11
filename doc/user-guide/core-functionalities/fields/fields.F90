program main

use atlas_module
implicit none

integer, parameter            :: wp = 8
integer                       :: jnode
character(len=:), allocatable :: unitsW, unitsP
type(atlas_Field)             :: pressureField , windField
type(atlas_Field)             :: pressureField2, windField2
real(wp), pointer             :: pressure(:), wind    (:,:)
type(atlas_FieldSet)          :: fields
type(atlas_metadata)          :: metadata
call atlas_init()

! Define fields
pressureField = atlas_Field("pressure", atlas_real(wp), [100])
windField     = atlas_Field("wind"    , atlas_real(wp), [2, 100])

! Initialize fields
call pressureField%data(pressure)
call windField    %data(wind)

! Assign values to fields
do jnode=1,100
  pressure(jnode) = 101325._wp
  wind(1, jnode)  = 0.01_wp + real(jnode,kind=wp)
  wind(2, jnode)  = 0.02_wp + real(jnode,kind=wp)
enddo

! Add info to fields
metadata = pressureField%metadata()
call metadata%set("units", "[Pa]")
call metadata%get("units", unitsP)
metadata = windField%metadata()
call metadata%set("units", "[m/s]")
call metadata%get("units", unitsW)

! Define fieldSet
fields = atlas_FieldSet("")
call fields%add(pressureField) ! Add pressureField to fieldSet
call fields%add(windField)     ! Add windField to fieldSet

! Retrieve field from fieldSet
pressureField2 = fields%field("pressure")
windField2     = fields%field("wind")

! Print some useful info
write(6, *) "name   = ", windField%name()
write(6, *) "size   = ", windField%size()
write(6, *) "units  = ", unitsW
write(6, *) "rank   = ", windField%rank()
write(6, *) "shape  = ", windField%shape(1), windField%shape(2)
write(6, *) "shape  = ", windField%shape()
write(6, *) "memory = ", windField%bytes(), "bytes"
write(6, *) "type   = ", windField%datatype()
write(6, *) "kind   = ", windField%kind()

! Print some values
write(6, *) "pressure(10) = ", pressure(10)
write(6, *) "wind(1, 10)  = ", wind(1,10)
write(6, *) "wind(2, 10)  = ", wind(2,10)

! Destroy object and free memory
call pressureField%final()
call windField    %final()
call fields       %final()

call atlas_finalize()
end program main
