!!$ This module sets up the floating point precision used.

module real_type_mod

!!$ sp contains the kind value for single precision.
  integer, parameter :: sp = kind(1.0)

!!$ dp contains the kind value for double precision.
  integer, parameter :: dp = kind(1.0d0)

!!$ Set wp to the desired precision.
  integer, parameter :: wp = dp

end module real_type_mod
