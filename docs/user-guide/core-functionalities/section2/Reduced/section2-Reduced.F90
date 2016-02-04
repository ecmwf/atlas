program run

use atlas_module

character(len=1024)     :: gridID
type(atlas_ReducedGrid) :: reducedGrid
real                    :: crdLat, crdLon
real                    :: crdLonLat(2)

call atlas_init()

call atlas_resource("--grid", "N32", gridID)

reducedGrid = atlas_ReducedGrid(gridID)

write(unit=6, fmt = "(A,I5.2,1x)") "nlat = ", reducedGrid%nlat()
write(unit=6, fmt = "(A,I5.2,1x)") "nlon = ", reducedGrid%nlon(0)
write(unit=6, fmt = "(A,I8.2,1x)") "npts = ", reducedGrid%npts()

crdLat = reducedGrid%lat(0)
crdLon = reducedGrid%lon(0, 1)
!crdLonLat = reducedGrid%lonlat(0, 1);

write(unit=6, fmt = "(A,F8.4,1x)") "lat  = ", crdLat
write(unit=6, fmt = "(A,F8.4,1x)") "lon  = ", crdLon

call atlas_finalize()

end program run
