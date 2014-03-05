program main
use common_module
use datastruct
use ftnunit
use iso_c_binding
use parallel_module
implicit none
type(Mesh_type) :: mesh
type(FunctionSpace_type) :: func_space
type(Field_type) :: field

call runtests_init
call parallel_init()

mesh = new_Mesh()
call runtests( test_all )
call delete(mesh)
call parallel_finalise()
call runtests_final

contains

subroutine test_all
  call test( test_function_space, "FunctionSpace" )
  call test( test_field_name, "Field_name" )
  call test( test_field_size, "Field_size" )
  call test( test_create_remove, "Field_create_remove" )
  call test( test_field_metadata, "Field_metadata" )
  call test( test_prismatic_function_space, "Prismatic" )
  call test( test_fieldset, "FieldSet" )
end subroutine

subroutine test_function_space
  call mesh%add_function_space( new_FunctionSpace("nodes", "P1-C", 10) )
  func_space = mesh%function_space("nodes")
  call assert_equal( func_space%name() , "nodes", "check_name" )
end subroutine

subroutine test_field_name
  call func_space%create_field("field",1,real_kind(jprb))
  field = func_space%field("field")
  call assert_equal( field%name() , "field", "check_name" )
end subroutine

subroutine test_field_metadata
  integer :: int
  logical :: true, false
  real(c_float) :: real32
  real(c_double) :: real64
  character(len=:), allocatable :: string
  type(MetaData_type) metadata
  call func_space%create_field("field_prop",1,real_kind(jprs))
  field = func_space%field("field_prop")

  metadata = field%metadata()

  call metadata%add("true",.True.)
  call metadata%add("false",.False.)
  call metadata%add("int",20)
  call metadata%add("real32", real(0.1,kind=jprs)   )
  call metadata%add("real64", real(0.2,kind=jprb) )
  call metadata%add("string", "hello world")

  call metadata%get("true",true)
  call metadata%get("false",false)
  call metadata%get("int",int)
  call metadata%get("real32",real32)
  call metadata%get("real64",real64)
  call metadata%get("string",string)
  
  call assert_equal( true, .True., "check bool property")
  call assert_equal( false, .False., "check bool property")
  call assert_equal( int, 20, "check int property")
  call assert_comparable( real32, real(0.1,kind=jprs), real(0.,kind=jprs), "check real32 property")
  call assert_comparable( real64, real(0.2,kind=jprb), real(0.,kind=jprb), "check real64 property")
  call assert_equal( string, "hello world", "check string property")
end subroutine



subroutine test_field_size
  integer, pointer :: fdata(:)
  call func_space%create_field("field_0",0,integer_kind())
  field = func_space%field("field_0")
  call field%access_data(fdata)
  call assert_equal( field%data_type() , "int32" , "check_data_type_0" )
  call assert_equal( size(fdata) , 0 , "check_size_0" )

  call func_space%create_field("field_1",1,integer_kind())
  field = func_space%field("field_1")
  call field%access_data(fdata)
  call assert_equal( size(fdata) , 10 , "check_size_1" )

  call func_space%create_field("field_2",2,integer_kind())
  field = func_space%field("field_2")
  call field%access_data(fdata)
  call assert_equal( size(fdata) , 20 , "check_size_2" )

end subroutine

subroutine test_create_remove
  call func_space%create_field("bla",1,integer_kind())
  field = func_space%field("bla")
  call assert_equal( field%name(), "bla", "check_name" )

  call func_space%remove_field("bla")
end subroutine

subroutine test_prismatic_function_space
  real(kind=jprb), pointer :: scalar(:,:)
  real(kind=jprs), pointer :: vector(:,:,:)
  call mesh%add_function_space( new_PrismaticFunctionSpace("prismatic", "P1-C", 5, 10 ) ) 

  func_space = mesh%function_space("prismatic")
  call func_space%create_field("vector_field",3,real_kind(jprs))
  field = func_space%field("vector_field")
  call field%access_data(vector)
  call assert_equal( size(vector),   150, "check_size" )
  call assert_equal( size(vector,1), 3,   "check_shape1" )
  call assert_equal( size(vector,2), 5,  "check_shape2" )
  call assert_equal( size(vector,3), 10,   "check_shape3" )

  call func_space%create_field("scalar_field",1,real_kind(jprb))
  field = func_space%field("scalar_field")
  scalar => field%data2()
  call assert_equal( size(scalar),   50, "check_size" )
  call assert_equal( size(scalar,1), 5,  "check_shape1" )
  call assert_equal( size(scalar,2), 10, "check_shape2" )
end subroutine

subroutine test_fieldset
  type(FieldSet_type) :: fieldset
  type(Field_type), allocatable :: fields(:)
  fieldset = new_FieldSet("fieldset")
  func_space = mesh%function_space("nodes")

  call fieldset%add_field( func_space%field("field_0") )
  call fieldset%add_field( func_space%field("field_1") )
  call fieldset%add_field( func_space%field("field_2") )

  func_space = mesh%function_space("prismatic")
  call fieldset%add_field( func_space%field("vector_field") )

  call assert_equal( fieldset%size(), 4, "check number of fields")

  field = fieldset%field(1)
  call assert_equal( field%name(), "field_0", "check name field_0")
  field = fieldset%field(2)
  call assert_equal( field%name(), "field_1", "check name field_1")
  field = fieldset%field(3)
  call assert_equal( field%name(), "field_2", "check name field_2")
  field = fieldset%field(4)
  call assert_equal( field%name(), "vector_field", "check name vector_field")


  call fieldset%get_array(fields)
  call assert_equal( size(fields), 4, "check array number of fields" )
  call assert_equal( fields(1)%name(), "field_0" , "check array name field_0")
  call assert_equal( fields(2)%name(), "field_1" , "check array name field_1")
  call assert_equal( fields(3)%name(), "field_2" , "check array name field_2")
  call assert_equal( fields(4)%name(), "vector_field" , "check array name vector_field")

  call delete(fieldset)
end subroutine

end program main

