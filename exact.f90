!!$Set the exact solution

subroutine exact_solution(nv, np, ngz, u, x, t)

  use real_type_mod
  use grid_parameters_mod
  use system_parameters_mod
  use method_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(out) :: u
  real(kind=wp), dimension(1-ngz:np+ngz), intent(in) :: x
  real(kind=wp), intent(in) :: t

!!$  The renormalized coordinates

  real(kind=wp), dimension(1-ngz:np+ngz) :: xbar
  integer :: crossing_number
  real(kind=wp) :: crossing_fraction

!!$  Constants

  real(kind=wp) :: pi

!!$  Scalar wave functions

  real(kind=wp), dimension(:), allocatable :: q1, q2, dq1, dq2

  integer :: ierr

  pi = 4.0_wp * atan(1.0_wp)

!!$  Work out the coordinate correction

  crossing_number = t / crossing_time
  crossing_fraction = t / crossing_time - crossing_number
  xbar = x - crossing_fraction * (xmax - xmin)
  where (xbar > xmax)
    xbar = xbar - (xmax - xmin)
  elsewhere (xbar < xmin)
    xbar = xbar + (xmax - xmin)
  end where

!!$  Set initial data

  if (advection) then

    if (advection_sine) then

      u(1,:) = sin(pi * (xbar))

    else if (advection_squarewave) then

      where (abs(xbar) < 1.0_wp / 3.0_wp)
        u(1,:) =  1.0_wp
      elsewhere
        u(1,:) = 0.0_wp
      end where

    else if (advection_bell) then
      where (abs(xbar) < 0.5_wp)
        u(1,:) = 1.0_wp - 4.0_wp * xbar**2
      elsewhere
        u(1,:) = 0.0_wp
      end where
  
    end if
    
  else if (burgers) then

!!$    The double shock tube. See thesis, page 13
!!$    Not certain that all the cases are correct yet.
    
    if (burgers_ql < burgers_qc) then

      if (burgers_qc < burgers_qr) then

        where (x < burgers_xl + burgers_ql * t)
          u(1,:) = burgers_ql
        elsewhere (x < burgers_xl + burgers_qc * t)
          u(1,:) = (x - burgers_xl) / t
        elsewhere (x < burgers_xr + burgers_qc * t)
          u(1,:) = burgers_qc
        elsewhere
          u(1,:) = burgers_qr
        end where

      else

        if (burgers_xl + burgers_qc * t > &
            0.5_wp * (burgers_qc + burgers_qr) * t) then
          write(*,*) 'Cannot use this data at time ', t
        end if
        
        where (x < burgers_xl + burgers_ql * t)
          u(1,:) = burgers_ql
        elsewhere (x < burgers_xl + burgers_qc * t)
          u(1,:) = (x - burgers_xl) / t
        elsewhere (x < burgers_xr + 0.5_wp * (burgers_qc + burgers_qr) * t)
          u(1,:) = burgers_qc
        elsewhere
          u(1,:) = burgers_qr
        end where

      end if
      
    else

      if (burgers_qc < burgers_qr) then

        where (x < burgers_xl + 0.5_wp * (burgers_ql + burgers_qc) * t)
          u(1,:) = burgers_ql
        elsewhere (x < burgers_xr + burgers_qc * t)
          u(1,:) = burgers_qc
        elsewhere (x < burgers_xr + burgers_qr * t)
          u(1,:) = (x - burgers_xr) / t
        elsewhere
          u(1,:) = burgers_qr
        end where

      else
        
        if (burgers_xl + 0.5_wp * (burgers_ql + burgers_qc) * t > &
            burgers_xr + 0.5_wp * (burgers_qc + burgers_qr) * t) then
          write(*,*) 'Cannot use this data at time ', t
          stop
        end if
        
        where (x < burgers_xl + 0.5_wp * (burgers_ql + burgers_qc) * t)
          u(1,:) = burgers_ql
        elsewhere (x < burgers_xr + 0.5_wp * (burgers_qc + burgers_qr) * t)
          u(1,:) = burgers_qc
        elsewhere
          u(1,:) = burgers_qr
        end where

      end if
      
    end if

  else if (scalarwave) then
    
    allocate(q1 (1-ngz:np+ngz),&
             q2 (1-ngz:np+ngz),&
             dq1(1-ngz:np+ngz),&
             dq2(1-ngz:np+ngz),&
             STAT = ierr)

    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate arrays!'
      STOP
    end if

    if (scalarwave_tsgaussian) then
    
      q1 = 0.5_wp * exp( -( (x - t) / scalarwave_sigma )**2 )
      q2 = 0.5_wp * exp( -( (x + t) / scalarwave_sigma )**2 )
      
      dq1 = -(x - t) / scalarwave_sigma**2 * &
           exp( -( (x - t) / scalarwave_sigma )**2 )
      dq2 = -(x + t) / scalarwave_sigma**2 * &
           exp( -( (x + t) / scalarwave_sigma )**2 )

    else if (scalarwave_separatinggaussian) then

      where ((x-t) > 0.0_wp)
        
        q1 = 0.5_wp * exp( -( (x - t) / scalarwave_sigma )**2 )
        dq1 = -(x - t) / scalarwave_sigma**2 * &
             exp( -( (x - t) / scalarwave_sigma )**2 )
        
      elsewhere
        
        q1  = 0.5_wp
        dq1 = 0.0_wp

      end where

      where ((x+t) < 0.0_wp)
    
        q2 = 0.5_wp * exp( -( (x + t) / scalarwave_sigma )**2 )
        dq2 = -(x + t) / scalarwave_sigma**2 * &
             exp( -( (x + t) / scalarwave_sigma )**2 )

      elsewhere
      
        q2  = 0.5_wp 
        dq2 = 0.0_wp

      end where

    end if
    
    u(1, :) =   q1 + q2
    u(2, :) = -dq1 + dq2
    u(3, :) =  dq1 + dq2

    deallocate(q1, q2, dq1, dq2)

  else if (euler) then

!!$    Possible two barrier shock problem

    where (x < euler_x_l)

      u(4, :) = euler_rho_l
      u(5, :) = euler_v_l
      u(6, :) = euler_p_l

    elsewhere (x < euler_x_r)

      u(4, :) = euler_rho_m
      u(5, :) = euler_v_m
      u(6, :) = euler_p_m

    elsewhere

      u(4, :) = euler_rho_r
      u(5, :) = euler_v_r
      u(6, :) = euler_p_r

    end where
    
    call primitive2conservative(nv, np, ngz, u)

  end if ! if advection/burgers/scalarwave/euler
  
end subroutine exact_solution
