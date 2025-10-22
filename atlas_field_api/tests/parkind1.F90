module parkind1
#ifdef SINGLE_PRECISION
integer, parameter :: JPRB = 4
#else
integer, parameter :: JPRB = 8
#endif
end module
