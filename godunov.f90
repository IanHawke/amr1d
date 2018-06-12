!!$Driver program for the Godunov solver

program godunov

  use real_type_mod
!!$  use parameters_mod
  use grid_parameters_mod
  use system_parameters_mod
  use method_parameters_mod
  
  implicit none

!!$  The grid functions; storage and pointers

  real(kind=wp), allocatable, dimension(:), target :: u1, u2
  real(kind=wp), dimension(:), pointer :: u, u_p, u_tmp

!!$  The flux function

  real(kind=wp), allocatable, dimension(:) :: flux

!!$  The spatial coordinate

  real(kind=wp), allocatable, dimension(:) :: x

!!$  The exact data

  real(kind=wp), allocatable, dimension(:) :: u_exact

!!$  The error

  real(kind=wp), allocatable, dimension(:) :: error

!!$  The current time

  real(kind=wp) :: t

!!$  Indices

  integer :: i, n, nconv

!!$  Error check

  integer :: ierr

!!$  Read parameters

  call read_parameters

  do nconv = 1, n_convergence

    n_points = n_points * 2
    dx = dx / 2.0_wp
    dt = dt / 2.0_wp
  
!!$  Now we actually do things

    t = 0.0_wp

    allocate(u1     (1-ghostzones:n_points+ghostzones),&
             u2     (1-ghostzones:n_points+ghostzones),&
             flux   (1-ghostzones:n_points+ghostzones),&
             x      (1-ghostzones:n_points+ghostzones),&
             u_exact(1-ghostzones:n_points+ghostzones),&
             error  (1-ghostzones:n_points+ghostzones),&
             STAT=ierr)

    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate arrays!'
      STOP
    end if
  
    do i = 1-ghostzones, n_points+ghostzones
      x(i) = xmin + (dble(i)-0.5_wp) * dx
    end do
    
    call exact_solution(u_exact, x, t)
    
    u1 = u_exact
    u2 = u_exact
    
    u   => u1
    u_p => u2

!!$  Having initialized, run until the end time 

    do while ( last_time - t > 0.000001_wp * dt)

      u_tmp => u
      u     => u_p
      u_p   => u_tmp
      
      t = t + dt

      call flux_formula(u_p, u_p, flux)

      do i = 1, n_points
        
        u(i) = u_p(i) + dtodx * (flux(i) - flux(i+1))
        
      end do
      
      call boundary_conditions(u)

    end do
      
    call exact_solution(u_exact, x, t)

    error = abs(u - u_exact)

    call output('godunov',x, u, u_exact, error)

    call output_errors('godunov',nconv, error)

    deallocate(u1     ,&
               u2     ,&
               flux   ,&
               x      ,&
               u_exact,&
               error    )

  end do
  
end program godunov
