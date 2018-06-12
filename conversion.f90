!!$Routines to convert between conserved and primitive variables

subroutine primitive2conservative(nv, np, ngz, u)

  use real_type_mod
  use system_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(inout) :: u

  if (euler) then

    u(1, :) = u(4, :)
    u(2, :) = u(4, :) * u(5, :)
    u(3, :) = u(6, :) / (euler_gamma - 1.0_wp) + &
         0.5_wp * u(4, :) * u(5, :)**2

  end if
  
end subroutine primitive2conservative

subroutine conservative2primitive(nv, np, ngz, u)

  use real_type_mod
  use system_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(inout) :: u

  integer :: i

  if (euler) then

    u(4, :) = u(1, :)
    u(5, :) = u(2, :) / u(1, :)
    u(6, :) =  (euler_gamma - 1.0_wp) * &
         (u(3, :) - 0.5_wp * u(2, :)**2 / u(1, :))

!!$    A hacked in floor

!!$    do i = 1 - ngz, np + ngz
!!$
!!$      if (u(4, i) < 0.0_wp) then
!!$        u(4, i) = 1.d-20
!!$        write(*,*) 'Floor on rho'
!!$        u(5, i) = u(2, i) / u(1, i)
!!$        u(6, i) =  (euler_gamma - 1.0_wp) * &
!!$             (u(3, i) - 0.5_wp * u(2, i)**2 / u(1, i))
!!$      end if
!!$      if (u(6, i) < 0.0_wp) then
!!$        u(6, i) = 1.d-20
!!$        write(*,*) 'Floor on p'
!!$      end if
!!$    
!!$    end do
    
  end if

end subroutine conservative2primitive

