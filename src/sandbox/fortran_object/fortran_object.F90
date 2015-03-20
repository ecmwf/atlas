! Purpose of this program is to test if the compiler can support
! the final keyword for specifying a destructor for a derived type.
! Unfortunately gfortran version < 4.9 does not support it,
!
! This means a destructor needs to be called manually,
! and be public.

module atlas_sandbox_fortran_object_module
   implicit none
   private
   public :: AnimalType, AnimalType__dtor

   type :: AnimalType
      private
      integer :: m_age
   contains
      procedure :: age
      procedure :: speak

! Finalization may not be supported by compiler
!    Only since gfortran >= 4.9

      !final :: AnimalType__dtor

   end type

! Declare constructor as interface with same name as type
   interface AnimalType
      module procedure AnimalType__ctor
   end interface

contains

   function AnimalType__ctor(age) result(self)
      type(AnimalType) :: self
      integer :: age
      write(0,'(A)') "Constructor Animal"
      self%m_age = age
   end function

   subroutine AnimalType__dtor(self)
      type(AnimalType), intent(inout) :: self
      write(0,'(A)') "Destroying animal"
   end subroutine

   function age(self)
      class(AnimalType), intent(inout) :: self
      integer :: age
      age = self%m_age
   end function

   subroutine speak(self)
      class(AnimalType), intent(in) :: self

      write(0,'(A)') "Animal::speak not overridden"
   end subroutine

end

program atlas_sandbox_fortran_object
use atlas_sandbox_fortran_object_module
implicit none

  type(AnimalType) :: animal

  animal = AnimalType(8)


  write(0,'(A,I0)') "age = ",animal%age()

  call animal%speak()

  ! In case Finalization is supported by compiler, following
  ! statement is not required
  call AnimalType__dtor(animal)

end program
