
program example_fortran 
use atlas_module
implicit none


call atlas_log%set_debug(10)


call atlas_log%info("Atlas initialized")
call atlas_log%info("version = ["//atlas_version()//"]")

call atlas_init()

call atlas_log%disconnect_stdout()
call atlas_log%connect_fortran_unit(6)

call atlas_log%channel_info %set_prefix("[%P] info  -- ")
call atlas_log%channel_debug%set_prefix("[%P] debug -- ")
call atlas_log%channel_stats%set_prefix("(S) -- ")

call atlas_log%debug("Here is some debugging information")
call atlas_log%stats("value = ...")

call atlas_log%info("Section")
call atlas_log%indent("|   ")
call atlas_log%info("Subsection")
call atlas_log%indent("|   ")
call atlas_log%info("Subsubsection")
call atlas_log%indent("|   ")
call atlas_log%info("stuff")
call atlas_log%debug("more stuff")
call atlas_log%stats("lolo")
call atlas_log%dedent()
call atlas_log%info("Subsubsection")
call atlas_log%dedent()
call atlas_log%info("Subsection")
call atlas_log%dedent()
call atlas_log%info("Section")
!call atlas_finalize()
call atlas_log%debug("Atlas finalized",lvl=0,flush=.True.)
write(6,'(A)') "exit"
end program
