! typy.f90
! Specification of real / integer kinds and character kinds.

module typy

  integer, parameter, public :: lowikind = selected_int_kind(5)

  integer, parameter, public :: dprec  = selected_real_kind(15,99)
  integer, parameter, public :: qprec  = selected_real_kind(30,999)
  integer, parameter, public :: sprec  = selected_real_kind(8,9)

  ! real / integer kinds
  integer, parameter, public :: rkind = selected_real_kind(15,99)
  integer, parameter, public :: ikind = selected_int_kind(10)
  integer, parameter, public :: likind = selected_int_kind(16)

  integer, parameter, public :: chkind = selected_char_kind("default")

end module typy
