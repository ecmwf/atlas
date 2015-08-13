! (C) Copyright 2013-2014 ECMWF.


function atlas_NodesFunctionSpace__mesh_halo(mesh,halo) result(function_space)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesFunctionSpace) :: function_space
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in), optional :: halo
  integer :: opt_halo
  opt_halo = 0
  if( present(halo) ) opt_halo = halo
  function_space%cpp_object_ptr = atlas__NodesFunctionSpace__new(c_str(""),mesh%cpp_object_ptr,opt_halo)
end function

function atlas_NodesFunctionSpace__name_mesh_halo(name,mesh,halo) result(function_space)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesFunctionSpace) :: function_space
  character(len=*), intent(in) :: name
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in), optional :: halo
  integer :: opt_halo
  opt_halo = 0
  if( present(halo) ) opt_halo = halo
  function_space%cpp_object_ptr = atlas__NodesFunctionSpace__new(c_str(name),mesh%cpp_object_ptr,opt_halo)
end function

subroutine atlas_NodesFunctionSpace__delete(this)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesFunctionSpace), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__NodesFunctionSpace__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine atlas_NodesFunctionSpace__delete

function atlas_NodesFunctionSpace__nb_nodes(this) result (nb_nodes)
  use atlas_nodesfunctionspace_c_binding
  integer :: nb_nodes
  class(atlas_NodesFunctionSpace), intent(in) :: this
  nb_nodes = atlas__NodesFunctionSpace__nb_nodes(this%cpp_object_ptr)
end function

function atlas_NodesFunctionSpace__mesh(this) result (mesh)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Mesh) :: mesh
  class(atlas_NodesFunctionSpace), intent(in) :: this
  mesh%cpp_object_ptr = atlas__NodesFunctionSpace__mesh(this%cpp_object_ptr)
end function

function atlas_NodesFunctionSpace__nodes(this) result (nodes)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Nodes) :: nodes
  class(atlas_NodesFunctionSpace), intent(in) :: this
  nodes%cpp_object_ptr = atlas__NodesFunctionSpace__nodes(this%cpp_object_ptr)
end function

function atlas_NodesFunctionSpace__create_field_kind(this,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  integer, intent(in) :: kind
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_field(this%cpp_object_ptr,c_str(""),kind)
end function

function atlas_NodesFunctionSpace__create_field_name_kind(this,name,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_field(this%cpp_object_ptr,c_str(name),kind)
end function

function atlas_NodesFunctionSpace__create_field_vars_kind(this,vars,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_field_vars(this%cpp_object_ptr,c_str(""),vars,size(vars),fortran_ordering,kind)
end function

function atlas_NodesFunctionSpace__create_field_name_vars_kind(this,name,vars,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_field_vars(this%cpp_object_ptr,c_str(name),vars,size(vars),fortran_ordering,kind)
end function

function atlas_NodesFunctionSpace__create_field_template(this,template) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: template
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_field_template(this%cpp_object_ptr,c_str(""),template%cpp_object_ptr)
end function

function atlas_NodesFunctionSpace__create_field_name_template(this,name,template) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: template
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_field_template(this%cpp_object_ptr,c_str(name),template%cpp_object_ptr)
end function



function atlas_NodesFunctionSpace__create_glb_field_kind(this,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  integer, intent(in) :: kind
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_global_field(this%cpp_object_ptr,c_str(""),kind)
end function

function atlas_NodesFunctionSpace__create_glb_field_name_kind(this,name,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_global_field(this%cpp_object_ptr,c_str(name),kind)
end function

function atlas_NodesFunctionSpace__create_glb_field_vars_kind(this,vars,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_global_field_vars(this%cpp_object_ptr,c_str(""),vars,size(vars),fortran_ordering,kind)
end function

function atlas_NodesFunctionSpace__create_glb_field_name_vars_kind(this,name,vars,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_global_field_vars(this%cpp_object_ptr,c_str(name),vars,size(vars),fortran_ordering,kind)
end function

function atlas_NodesFunctionSpace__create_glb_field_template(this,template) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: template
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_global_field_template(this%cpp_object_ptr,c_str(""),template%cpp_object_ptr)
end function

function atlas_NodesFunctionSpace__create_glb_field_name_template(this,name,template) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: template
  field%cpp_object_ptr = atlas__NodesFunctionSpace__create_global_field_template(this%cpp_object_ptr,c_str(name),template%cpp_object_ptr)
end function

subroutine atlas_NodesFunctionSpace__halo_exchange_fieldset(this,fieldset)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_FieldSet), intent(inout) :: fieldset
  call atlas__NodesFunctionSpace__halo_exchange_fieldset(this%cpp_object_ptr,fieldset%cpp_object_ptr)
end subroutine

subroutine atlas_NodesFunctionSpace__halo_exchange_field(this,field)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field), intent(inout) :: field
  call atlas__NodesFunctionSpace__halo_exchange_field(this%cpp_object_ptr,field%cpp_object_ptr)
end subroutine

subroutine atlas_NodesFunctionSpace__gather_fieldset(this,local,global)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: local
  type(atlas_FieldSet), intent(inout) :: global
  call atlas__NodesFunctionSpace__gather_fieldset(this%cpp_object_ptr,local%cpp_object_ptr,global%cpp_object_ptr)
end subroutine

subroutine atlas_NodesFunctionSpace__gather_field(this,local,global)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__NodesFunctionSpace__gather_field(this%cpp_object_ptr,local%cpp_object_ptr,global%cpp_object_ptr)
end subroutine

subroutine atlas_NodesFunctionSpace__scatter_fieldset(this,global,local)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: global
  type(atlas_FieldSet), intent(inout) :: local
  call atlas__NodesFunctionSpace__scatter_fieldset(this%cpp_object_ptr,global%cpp_object_ptr,local%cpp_object_ptr)
end subroutine

subroutine atlas_NodesFunctionSpace__scatter_field(this,global,local)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__NodesFunctionSpace__scatter_field(this%cpp_object_ptr,global%cpp_object_ptr,local%cpp_object_ptr)
end subroutine

function atlas_NodesFunctionSpace__checksum_fieldset(this,fieldset) result(checksum)
  use atlas_nodesfunctionspace_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__NodesFunctionSpace__checksum_fieldset(this%cpp_object_ptr,fieldset%cpp_object_ptr,checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function


function atlas_NodesFunctionSpace__checksum_field(this,field) result(checksum)
  use atlas_nodesfunctionspace_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__NodesFunctionSpace__checksum_field(this%cpp_object_ptr,field%cpp_object_ptr,checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function


subroutine atlas_NodesFunctionSpace__minimum_real32_r0(this,field,minimum)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: minimum
  call atlas__NodesFunctionSpace__min_float(this%cpp_object_ptr,field%cpp_object_ptr,minimum)
end subroutine

subroutine atlas_NodesFunctionSpace__minimum_real32_r1(this,field,minimum)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: minimum(:)
  type(c_ptr) :: min_cptr
  real(c_float), pointer :: min_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__min_arr_float(this%cpp_object_ptr,field%cpp_object_ptr,min_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  allocate(minimum(min_size))
  minimum(:) = min_fptr(:)
  call atlas_free(min_cptr)
end subroutine

subroutine atlas_NodesFunctionSpace__maximum_real32_r0(this,field,maximum)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: maximum
  call atlas__NodesFunctionSpace__max_float(this%cpp_object_ptr,field%cpp_object_ptr,maximum)
end subroutine

subroutine atlas_NodesFunctionSpace__maximum_real32_r1(this,field,maximum)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: maximum(:)
  type(c_ptr) :: max_cptr
  real(c_float), pointer :: max_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__max_arr_float(this%cpp_object_ptr,field%cpp_object_ptr,max_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  allocate(maximum(max_size))
  maximum(:) = max_fptr(:)
  call atlas_free(max_cptr)
end subroutine


subroutine atlas_NodesFunctionSpace__minloc_real32_r0(this,field,minimum,location)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__minloc_float(this%cpp_object_ptr,field%cpp_object_ptr,minimum,loc)
  location = loc
end subroutine

subroutine atlas_NodesFunctionSpace__maxloc_real32_r0(this,field,maximum,location)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__maxloc_float(this%cpp_object_ptr,field%cpp_object_ptr,maximum,loc)
  location = loc
end subroutine


subroutine atlas_NodesFunctionSpace__minloc_real32_r1(this,field,minimum,location)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: min_cptr, loc_cptr
  real(c_float), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloc_arr_float(this%cpp_object_ptr,field%cpp_object_ptr,min_cptr,loc_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  call atlas_free(min_cptr)
  call atlas_free(loc_cptr)
end subroutine

subroutine atlas_NodesFunctionSpace__maxloc_real32_r1(this,field,maximum,location)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: max_cptr, loc_cptr
  real(c_float), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloc_arr_float(this%cpp_object_ptr,field%cpp_object_ptr,max_cptr,loc_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  call atlas_free(max_cptr)
  call atlas_free(loc_cptr)
end subroutine


subroutine atlas_NodesFunctionSpace__sum_real32_r0(this,field,sum,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_float(this%cpp_object_ptr,field%cpp_object_ptr,sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__sum_real32_r1(this,field,sum,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  real(c_float), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_arr_float(this%cpp_object_ptr,field%cpp_object_ptr,sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call atlas_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__order_independent_sum_real32_r0(this,field,sum,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_float(this%cpp_object_ptr,field%cpp_object_ptr,sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__order_independent_sum_real32_r1(this,field,sum,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  real(c_float), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_arr_float(this%cpp_object_ptr,field%cpp_object_ptr,sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call atlas_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__mean_real32_r0(this,field,mean,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: mean
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_float(this%cpp_object_ptr,field%cpp_object_ptr,mean,opt_N)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__mean_real32_r1(this,field,mean,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: mean(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr
  real(c_float), pointer :: mean_fptr(:)
  integer :: mean_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_arr_float(this%cpp_object_ptr,field%cpp_object_ptr,mean_cptr,mean_size,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/mean_size/))
  allocate(mean(mean_size))
  mean(:) = mean_fptr(:)
  call atlas_free(mean_cptr)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__mean_and_stddev_real32_r0(this,field,mean,stddev,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: mean
  real(c_float), intent(out) :: stddev
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_float(this%cpp_object_ptr,field%cpp_object_ptr,mean,stddev,opt_N)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__mean_and_stddev_real32_r1(this,field,mean,stddev,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: mean(:)
  real(c_float), allocatable, intent(out) :: stddev(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr, stddev_cptr
  real(c_float), pointer :: mean_fptr(:), stddev_fptr(:)
  integer :: varsize
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_arr_float(this%cpp_object_ptr,field%cpp_object_ptr,mean_cptr,stddev_cptr,varsize,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/varsize/))
  call c_f_pointer(stddev_cptr,stddev_fptr,(/varsize/))
  allocate(mean(varsize))
  allocate(stddev(varsize))
  mean(:) = mean_fptr(:)
  stddev(:) = stddev_fptr(:)
  call atlas_free(mean_cptr)
  call atlas_free(stddev_cptr)
  if( present(N) ) N = opt_N
end subroutine



subroutine atlas_NodesFunctionSpace__minimum_real64_r0(this,field,minimum)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: minimum
  call atlas__NodesFunctionSpace__min_double(this%cpp_object_ptr,field%cpp_object_ptr,minimum)
end subroutine

subroutine atlas_NodesFunctionSpace__minimum_real64_r1(this,field,minimum)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: minimum(:)
  type(c_ptr) :: min_cptr
  real(c_double), pointer :: min_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__min_arr_double(this%cpp_object_ptr,field%cpp_object_ptr,min_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  allocate(minimum(min_size))
  minimum(:) = min_fptr(:)
  call atlas_free(min_cptr)
end subroutine

subroutine atlas_NodesFunctionSpace__maximum_real64_r0(this,field,maximum)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: maximum
  call atlas__NodesFunctionSpace__max_double(this%cpp_object_ptr,field%cpp_object_ptr,maximum)
end subroutine

subroutine atlas_NodesFunctionSpace__maximum_real64_r1(this,field,maximum)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: maximum(:)
  type(c_ptr) :: max_cptr
  real(c_double), pointer :: max_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__max_arr_double(this%cpp_object_ptr,field%cpp_object_ptr,max_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  allocate(maximum(max_size))
  maximum(:) = max_fptr(:)
  call atlas_free(max_cptr)
end subroutine


subroutine atlas_NodesFunctionSpace__minloc_real64_r0(this,field,minimum,location)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__minloc_double(this%cpp_object_ptr,field%cpp_object_ptr,minimum,loc)
  location = loc
end subroutine

subroutine atlas_NodesFunctionSpace__maxloc_real64_r0(this,field,maximum,location)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__maxloc_double(this%cpp_object_ptr,field%cpp_object_ptr,maximum,loc)
  location = loc
end subroutine


subroutine atlas_NodesFunctionSpace__minloc_real64_r1(this,field,minimum,location)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: min_cptr, loc_cptr
  real(c_double), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloc_arr_double(this%cpp_object_ptr,field%cpp_object_ptr,min_cptr,loc_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  call atlas_free(min_cptr)
  call atlas_free(loc_cptr)
end subroutine

subroutine atlas_NodesFunctionSpace__maxloc_real64_r1(this,field,maximum,location)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: max_cptr, loc_cptr
  real(c_double), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloc_arr_double(this%cpp_object_ptr,field%cpp_object_ptr,max_cptr,loc_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  call atlas_free(max_cptr)
  call atlas_free(loc_cptr)
end subroutine


subroutine atlas_NodesFunctionSpace__sum_real64_r0(this,field,sum,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_double(this%cpp_object_ptr,field%cpp_object_ptr,sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__sum_real64_r1(this,field,sum,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  real(c_double), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_arr_double(this%cpp_object_ptr,field%cpp_object_ptr,sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call atlas_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__order_independent_sum_real64_r0(this,field,sum,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_double(this%cpp_object_ptr,field%cpp_object_ptr,sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__order_independent_sum_real64_r1(this,field,sum,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  real(c_double), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_arr_double(this%cpp_object_ptr,field%cpp_object_ptr,sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call atlas_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__mean_real64_r0(this,field,mean,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: mean
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_double(this%cpp_object_ptr,field%cpp_object_ptr,mean,opt_N)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__mean_real64_r1(this,field,mean,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: mean(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr
  real(c_double), pointer :: mean_fptr(:)
  integer :: mean_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_arr_double(this%cpp_object_ptr,field%cpp_object_ptr,mean_cptr,mean_size,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/mean_size/))
  allocate(mean(mean_size))
  mean(:) = mean_fptr(:)
  call atlas_free(mean_cptr)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__mean_and_stddev_real64_r0(this,field,mean,stddev,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: mean
  real(c_double), intent(out) :: stddev
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_double(this%cpp_object_ptr,field%cpp_object_ptr,mean,stddev,opt_N)
  if( present(N) ) N = opt_N
end subroutine

subroutine atlas_NodesFunctionSpace__mean_and_stddev_real64_r1(this,field,mean,stddev,N)
use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesFunctionSpace), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: mean(:)
  real(c_double), allocatable, intent(out) :: stddev(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr, stddev_cptr
  real(c_double), pointer :: mean_fptr(:), stddev_fptr(:)
  integer :: varsize
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_arr_double(this%cpp_object_ptr,field%cpp_object_ptr,mean_cptr,stddev_cptr,varsize,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/varsize/))
  call c_f_pointer(stddev_cptr,stddev_fptr,(/varsize/))
  allocate(mean(varsize))
  allocate(stddev(varsize))
  mean(:) = mean_fptr(:)
  stddev(:) = stddev_fptr(:)
  call atlas_free(mean_cptr)
  call atlas_free(stddev_cptr)
  if( present(N) ) N = opt_N
end subroutine

function atlas_NodesColumnFunctionSpace__mesh_levels_halo(mesh,nb_levels,halo) result(function_space)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesColumnFunctionSpace) :: function_space
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in) :: nb_levels
  integer, intent(in), optional :: halo
  integer :: opt_halo
  opt_halo = 0
  if( present(halo) ) opt_halo = halo
  function_space%cpp_object_ptr = atlas__NodesColumnFunctionSpace__new(c_str(""),mesh%cpp_object_ptr,nb_levels,opt_halo)
end function

function atlas_NodesColumnFunctionSpace__name_mesh_levels_halo(name,mesh,nb_levels,halo) result(function_space)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesColumnFunctionSpace) :: function_space
  character(len=*), intent(in) :: name
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in) :: nb_levels
  integer, intent(in), optional :: halo
  integer :: opt_halo
  opt_halo = 0
  if( present(halo) ) opt_halo = halo
  function_space%cpp_object_ptr = atlas__NodesColumnFunctionSpace__new(c_str(name),mesh%cpp_object_ptr,nb_levels,opt_halo)
end function

subroutine atlas_NodesColumnFunctionSpace__delete(this)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesColumnFunctionSpace), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__NodesColumnFunctionSpace__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine

function atlas_NodesColumnFunctionSpace__nb_nodes(this) result (nb_nodes)
  use atlas_nodesfunctionspace_c_binding
  integer :: nb_nodes
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  nb_nodes = atlas__NodesFunctionSpace__nb_nodes(this%cpp_object_ptr)
end function

function atlas_NodesColumnFunctionSpace__nb_levels(this) result (nb_levels)
  use atlas_nodesfunctionspace_c_binding
  integer :: nb_levels
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  nb_levels = atlas__NodesColumnFunctionSpace__nb_levels(this%cpp_object_ptr)
end function

function atlas_NodesColumnFunctionSpace__mesh(this) result (mesh)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Mesh) :: mesh
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  mesh%cpp_object_ptr = atlas__NodesFunctionSpace__mesh(this%cpp_object_ptr)
end function

function atlas_NodesColumnFunctionSpace__nodes(this) result (nodes)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Nodes) :: nodes
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  nodes%cpp_object_ptr = atlas__NodesFunctionSpace__nodes(this%cpp_object_ptr)
end function


function atlas_NCFunctionSpace__create_field_kind(this,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  integer, intent(in) :: kind
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_field(this%cpp_object_ptr,c_str(""),kind)
end function

function atlas_NCFunctionSpace__create_field_name_kind(this,name,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_field(this%cpp_object_ptr,c_str(name),kind)
end function

function atlas_NCFunctionSpace__create_field_vars_kind(this,vars,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_field_vars(this%cpp_object_ptr,c_str(""),vars,size(vars),fortran_ordering,kind)
end function

function atlas_NCFunctionSpace__create_field_name_vars_kind(this,name,vars,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_field_vars(this%cpp_object_ptr,c_str(name),vars,size(vars),fortran_ordering,kind)
end function

function atlas_NCFunctionSpace__create_field_template(this,template) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_Field) :: template
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_field_template(this%cpp_object_ptr,c_str(""),template%cpp_object_ptr)
end function

function atlas_NCFunctionSpace__create_field_name_template(this,name,template) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: template
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_field_template(this%cpp_object_ptr,c_str(name),template%cpp_object_ptr)
end function

function atlas_NCFunctionSpace__create_glb_field_kind(this,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  integer, intent(in) :: kind
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_global_field(this%cpp_object_ptr,c_str(""),kind)
end function

function atlas_NCFunctionSpace__create_glb_field_name_kind(this,name,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_global_field(this%cpp_object_ptr,c_str(name),kind)
end function

function atlas_NCFunctionSpace__create_glb_field_vars_kind(this,vars,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_global_field_vars(this%cpp_object_ptr,c_str(""),vars,size(vars),fortran_ordering,kind)
end function

function atlas_NCFunctionSpace__create_glb_field_name_vars_kind(this,name,vars,kind) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_global_field_vars(this%cpp_object_ptr,c_str(name),vars,size(vars),fortran_ordering,kind)
end function

function atlas_NCFunctionSpace__create_glb_field_template(this,template) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_Field) :: template
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_global_field_template(this%cpp_object_ptr,c_str(""),template%cpp_object_ptr)
end function

function atlas_NCFunctionSpace__create_glb_field_name_template(this,name,template) result(field)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: template
  field%cpp_object_ptr = atlas__NodesColumnFunctionSpace__create_global_field_template(this%cpp_object_ptr,c_str(name),template%cpp_object_ptr)
end function


subroutine atlas_NCFunctionSpace__halo_exchange_fieldset(this,fieldset)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_FieldSet), intent(inout) :: fieldset
  call atlas__NodesFunctionSpace__halo_exchange_fieldset(this%cpp_object_ptr,fieldset%cpp_object_ptr)
end subroutine

subroutine atlas_NCFunctionSpace__halo_exchange_field(this,field)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_Field), intent(inout) :: field
  call atlas__NodesFunctionSpace__halo_exchange_field(this%cpp_object_ptr,field%cpp_object_ptr)
end subroutine

subroutine atlas_NCFunctionSpace__gather_fieldset(this,local,global)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: local
  type(atlas_FieldSet), intent(inout) :: global
  call atlas__NodesFunctionSpace__gather_fieldset(this%cpp_object_ptr,local%cpp_object_ptr,global%cpp_object_ptr)
end subroutine

subroutine atlas_NCFunctionSpace__gather_field(this,local,global)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__NodesFunctionSpace__gather_field(this%cpp_object_ptr,local%cpp_object_ptr,global%cpp_object_ptr)
end subroutine

subroutine atlas_NCFunctionSpace__scatter_fieldset(this,global,local)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: global
  type(atlas_FieldSet), intent(inout) :: local
  call atlas__NodesFunctionSpace__scatter_fieldset(this%cpp_object_ptr,global%cpp_object_ptr,local%cpp_object_ptr)
end subroutine

subroutine atlas_NCFunctionSpace__scatter_field(this,global,local)
  use atlas_nodesfunctionspace_c_binding
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__NodesFunctionSpace__scatter_field(this%cpp_object_ptr,global%cpp_object_ptr,local%cpp_object_ptr)
end subroutine

function atlas_NCFunctionSpace__checksum_fieldset(this,fieldset) result(checksum)
  use atlas_nodesfunctionspace_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__NodesColumnFunctionSpace__checksum_fieldset(this%cpp_object_ptr,fieldset%cpp_object_ptr,checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function


function atlas_NCFunctionSpace__checksum_field(this,field) result(checksum)
  use atlas_nodesfunctionspace_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_NodesColumnFunctionSpace), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__NodesColumnFunctionSpace__checksum_field(this%cpp_object_ptr,field%cpp_object_ptr,checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function


