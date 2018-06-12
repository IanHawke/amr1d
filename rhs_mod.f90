!!$Compute the RHS function using TVD and exact flux

module rhs_mod

  use grid_mod
  use sources_mod

  contains

    subroutine compute_RHS(ngz, u, one_grid)
  
      implicit none

      integer, intent(in) :: ngz
      real(kind=wp), dimension(:, 1-ngz:), intent(in) :: u
      type(grid) :: one_grid

      real(kind=wp), dimension(:,:), pointer :: rhs, flux

      real(kind=wp), dimension(:,:), allocatable :: source

      integer :: ierr

      integer :: nv, np
      real(kind=wp) :: dx

      nv = one_grid%nvars
      np = one_grid%n_points
      dx = one_grid%dx

      allocate(source(nv, 1-ngz:np+ngz), &
               STAT = ierr)
      if (ierr .ne. 0) then
        write(*,*) 'Failed to allocate arrays!'
        STOP
      end if
  
      rhs => one_grid%grid_data%rhs
      flux => one_grid%grid_data%flux

      if (centred_differencing) then

        call pointwise_flux(nv, np, ngz, u, flux)

        call derivative(np, ngz, dx, flux, rhs)

!!$    Because of the ordering of the equation, the sign is wrong

        rhs = -rhs

      else

        call intercell_flux(nv, np, ngz, u, one_grid%grid_data)
    
        rhs(:, 1:np) = &
             (flux(:, 1:np) - flux(:, 2:np+1)) / dx

      end if

      call sources(ngz, u, source)

      rhs = rhs + source

      deallocate(source)

    end subroutine compute_rhs

!!$Straightforward centred differencing operators

    subroutine derivative(np, ngz, dx, u, du)

      use real_type_mod
      use method_parameters_mod

      implicit none

      integer, intent(in) :: np, ngz
      real(kind=wp), dimension(:, 1-ngz:), intent(in) :: u
      real(kind=wp), dimension(:, 1-ngz:), intent(out) :: du
      real(kind=wp), intent(in) :: dx

      integer :: i

      if (centred_differencing_order == 2) then

        do i = 1, np

          du(:, i) = (u(:, i+1) - u(:, i-1)) / (2.0_wp * dx)

        end do
    
      else if (centred_differencing_order == 4) then

        do i = 1, np

          du(:, i) = ( -u(:, i+2) + u(:, i-2) + &
               8.0_wp * (u(:, i+1) - u(:, i-1)) ) / (12.0_wp * dx)

        end do

      else

        write(*,*) 'Order of centred differencing not recognized!'
        stop

      end if
  
    end subroutine derivative

end module rhs_mod
