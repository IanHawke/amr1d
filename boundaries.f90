!!$Set the boundary conditions

subroutine boundary_conditions(nv, np, ngz, u)

  use real_type_mod
  use method_parameters_mod
  use system_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(inout) :: u

  integer :: i

  if (periodic) then

    u(:, 1-ngz:0)                 = u(:, np-ngz+1:np)
    u(:, np+1:np+ngz)             = u(:, 1:ngz)

  else if (outflow) then

    do i = 1, ngz

      u(:, 1 - i )                 = u(:, 1)
      u(:, np + i)                 = u(:, np)
      
    end do

  else if (reflection) then

    do i = 1, ngz

      u(:, 1 - i )                 = u(:, 1)
      u(:, np + i)                 = u(:, np)
      
    end do    
      
    if (euler) then

      do i = 1, ngz

        u(2, 1 - i )                 = -u(2, i)
        u(5, 1 - i )                 = -u(5, i)
        u(2, np + i)                 = -u(2, np + 1 - i)
        u(5, np + i)                 = -u(5, np + 1 - i)

      end do
      
    end if
  
  end if
  
end subroutine boundary_conditions

subroutine boundary_conditions_left(nv, np, ngz, u)

  use real_type_mod
  use method_parameters_mod
  use system_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(inout) :: u

  integer :: i

  if (periodic) then

    write(*,*) 'We cannot set just one boundary when using periodic!'
    STOP

  else if (outflow) then

    do i = 1, ngz

      u(:, 1 - i )                 = u(:, 1)
      
    end do

  else if (reflection) then

    do i = 1, ngz

      u(:, 1 - i )                 = u(:, 1)
      
    end do    
      
    if (euler) then

      do i = 1, ngz

        u(2, 1 - i )                 = -u(2, i)
        u(5, 1 - i )                 = -u(5, i)

      end do
      
    end if
  
  end if
  
end subroutine boundary_conditions_left

subroutine boundary_conditions_right(nv, np, ngz, u)

  use real_type_mod
  use method_parameters_mod
  use system_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(inout) :: u

  integer :: i

  if (periodic) then

    write(*,*) 'We cannot set just one boundary when using periodic!'
    STOP

  else if (outflow) then

    do i = 1, ngz

      u(:, np + i)                 = u(:, np)
      
    end do

  else if (reflection) then

    do i = 1, ngz

      u(:, np + i)                 = u(:, np)
      
    end do    
      
    if (euler) then

      do i = 1, ngz

        u(2, np + i)                 = -u(2, np + 1 - i)
        u(5, np + i)                 = -u(5, np + 1 - i)

      end do
      
    end if
  
  end if
  
end subroutine boundary_conditions_right
