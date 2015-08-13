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


