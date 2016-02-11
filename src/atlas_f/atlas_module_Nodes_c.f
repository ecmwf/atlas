! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
!

function Nodes__metadata(this) result(metadata)
  use atlas_nodes_c_binding
  type(atlas_Metadata) :: metadata
  class(atlas_mesh_Nodes), intent(in) :: this
  call metadata%reset_c_ptr( atlas__mesh__Nodes__metadata(this%c_ptr()) )
end function

function Nodes__size(this) result(size)
  use atlas_nodes_c_binding
  integer :: size
  class(atlas_mesh_Nodes), intent(in) :: this
  size = atlas__mesh__Nodes__size(this%c_ptr())
end function

subroutine Nodes__resize(this,size)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  integer, intent(in) :: size
  call atlas__mesh__Nodes__resize(this%c_ptr(),size)
end subroutine

subroutine Nodes__add(this,field)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(inout) :: this
  type(atlas_Field), intent(in) :: field
  call atlas__mesh__Nodes__add(this%c_ptr(),field%c_ptr())
end subroutine

subroutine Nodes__remove_field(this,name)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__mesh__Nodes__remove_field(this%c_ptr(),c_str(name))
end subroutine

function Nodes__has_field(this,name) result(has_field)
  use atlas_nodes_c_binding
  logical :: has_field
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  integer :: has_field_int
  has_field_int = atlas__mesh__Nodes__has_field(this%c_ptr(),c_str(name))
  has_field = .False.
  if( has_field_int == 1 ) has_field = .True.
end function

function Nodes__nb_fields(this) result(nb_fields)
  use atlas_nodes_c_binding
  integer :: nb_fields
  class(atlas_mesh_Nodes), intent(in) :: this
  nb_fields = atlas__mesh__Nodes__nb_fields(this%c_ptr())
end function

function Nodes__field_by_name(this,name) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function Nodes__field_by_idx(this,idx) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  integer, intent(in) :: idx
  field = atlas_Field( atlas__mesh__Nodes__field_by_idx(this%c_ptr(),idx) )
  call field%return()
end function

function Nodes__str(this) result(str)
  use atlas_nodes_c_binding
  character(len=:), allocatable :: str
  class(atlas_mesh_Nodes), intent(in) :: this
  type(c_ptr) :: str_cptr
  integer(c_int) :: str_size
  call atlas__mesh__Nodes__str(this%c_ptr(),str_cptr,str_size)
  allocate(character(len=str_size) :: str )
  str = c_to_f_string_cptr(str_cptr)
  call atlas_free(str_cptr)
end function


function Nodes__lonlat(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%c_ptr(),c_str("lonlat")) )
  call field%return()
end function

function Nodes__ghost(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%c_ptr(),c_str("ghost")) )
  call field%return()
end function

function Nodes__global_index(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%c_ptr(),c_str("glb_idx")) )
  call field%return()
end function

function Nodes__remote_index(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%c_ptr(),c_str("remote_idx")) )
  call field%return()
end function

function Nodes__partition(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%c_ptr(),c_str("partition")) )
  call field%return()
end function

subroutine atlas_mesh_Nodes__delete(this)
  class(atlas_mesh_Nodes), intent(inout) :: this
end subroutine

subroutine atlas_mesh_Nodes__copy(this,obj_in)
  class(atlas_mesh_Nodes), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine
