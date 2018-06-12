!!$Driver program for the slope limited TVD solver

program tvd

  use real_type_mod
  use grid_parameters_mod
  use grid_mod
!!$  use glist_mod
    
  implicit none

!!$  The grid functions; storage and pointers

!!$  real(kind=wp), dimension(:,:), pointer :: u, u_p
!!$  real(kind=wp), dimension(:),   pointer :: x

!!$  The spatial coordinate

!!$  real(kind=wp), allocatable, dimension(:) :: x

  type(grid), target :: base

!!$  The exact data

  real(kind=wp), allocatable, dimension(:,:) :: u_exact

!!$  The error

  real(kind=wp), allocatable, dimension(:,:) :: error

!!$  The current time

  real(kind=wp) :: t

!!$  Indices

  integer :: i, k, n, nconv

!!$  Error check

  integer :: ierr

!!$  Testing

  type(grid) :: child
  type(box)  :: child_box
  type(grid), pointer :: parent

!!$  type(glist) :: list

!!$  Read parameters

  call read_parameters

!!$  Check

  if (ghostzones < 2) then
    write(*,*) 'TVD requires 2 ghostzones!'
    stop
  end if

  do nconv = 1, n_convergence

    base_n_points = base_n_points * 2
    base_dx = base_dx / 2.0_wp
    base_dt = base_dt / 2.0_wp
  
!!$    Just testing

!!$    call setup_base_glist(list)
!!$    call delete_glist(list)

!!$    Testing prolongation and restriction for just one timelevel
!!$
!!$    t = 0.0_wp
!!$
!!$    call setup_base_grid(base)
!!$
!!$    allocate(u_exact(base%nvars, 1-ghostzones:base%n_points+ghostzones),&
!!$             error  (base%nvars, 1-ghostzones:base%n_points+ghostzones),&
!!$             STAT=ierr)
!!$
!!$    if (ierr .ne. 0) then
!!$      write(*,*) 'Failed to allocate arrays!'
!!$      STOP
!!$    end if
!!$
!!$    call exact_solution(base%nvars, base%n_points, ghostzones, u_exact, &
!!$         base%x, t)
!!$
!!$    base%u   = u_exact
!!$    base%u_p = u_exact

!!$    child_box%lbnd = base%mybox%lbnd
!!$    child_box%ubnd = 2 * base%mybox%ubnd
!!$    child_box%npoints = 2 * base%mybox%npoints
!!$    child_box%nghosts = base%mybox%nghosts
!!$    child_box%bbox = base%mybox%bbox
!!$    parent => base
!!$
!!$    child_box%lbnd = base%mybox%lbnd + 4
!!$    child_box%ubnd = base%mybox%ubnd - 4
!!$    child_box%npoints = 2 * (child_box%ubnd - child_box%lbnd + 1)
!!$    child_box%nghosts = base%mybox%nghosts
!!$    child_box%bbox = 0
!!$    parent => base
!!$
!!$    call setup_grid(child, child_box, &
!!$         base%nvars, &
!!$         base%dx / 2.0_wp, &
!!$         base%dt / 2.0_wp, &
!!$         parent)
!!$
!!$    call prolongate_grid(child)
!!$
!!$    do i = 1 - base%ghostzones, base%n_points + base%ghostzones
!!$      write(11, '(es21.12E3,a1)',advance='NO') base%x(i), ' '
!!$      do k = 1, base%nvars
!!$        write(11, '(es21.12E3,a1)',advance='NO') base%u(k,i), ' '
!!$      end do
!!$      write(11, *)
!!$    end do
!!$
!!$    do i = 1 - child%ghostzones, child%n_points + child%ghostzones
!!$      write(12, '(es21.12E3,a1)',advance='NO') child%x(i), ' '
!!$      do k = 1, child%nvars
!!$        write(12, '(es21.12E3,a1)',advance='NO') child%u(k,i), ' '
!!$      end do
!!$      write(12, *)
!!$    end do

!!$    where (abs(child%x) < 0.295_wp)
!!$
!!$      child%u(1, :) = 0.0_wp
!!$
!!$    end where
    
!!$    call restrict_grid(child)
!!$
!!$    do i = 1 - base%ghostzones, base%n_points + base%ghostzones
!!$      write(13, '(es21.12E3,a1)',advance='NO') base%x(i), ' '
!!$      do k = 1, base%nvars
!!$        write(13, '(es21.12E3,a1)',advance='NO') base%u(k,i), ' '
!!$      end do
!!$      write(13, *)
!!$    end do
!!$
!!$    do i = 1 - child%ghostzones, child%n_points + child%ghostzones
!!$      write(14, '(es21.12E3,a1)',advance='NO') child%x(i), ' '
!!$      do k = 1, child%nvars
!!$        write(14, '(es21.12E3,a1)',advance='NO') child%u(k,i), ' '
!!$      end do
!!$      write(14, *)
!!$    end do
!!$
!!$    write(*,*) 'Got this far'
!!$    stop

!!$    call swap_timelevels(base)
!!$    call evolve_one_step(base)
!!$
!!$    call swap_timelevels(child)
!!$    call evolve_one_step(child)
!!$    call swap_timelevels(child)
!!$    call evolve_one_step(child)
!!$
!!$    call restrict_grid(child)
!!$
!!$    call prolongate_boundaries(child, child%u, t)
!!$
!!$    do i = 1 - base%ghostzones, base%n_points + base%ghostzones
!!$      write(13, '(es21.12E3,a1)',advance='NO') base%x(i), ' '
!!$      do k = 1, base%nvars
!!$        write(13, '(es21.12E3,a1)',advance='NO') base%u(k,i), ' '
!!$      end do
!!$      write(13, *)
!!$    end do
!!$
!!$    do i = 1 - child%ghostzones, child%n_points + child%ghostzones
!!$      write(14, '(es21.12E3,a1)',advance='NO') child%x(i), ' '
!!$      do k = 1, child%nvars
!!$        write(14, '(es21.12E3,a1)',advance='NO') child%u(k,i), ' '
!!$      end do
!!$      write(14, *)
!!$    end do
!!$
!!$    write(*,*) 'Got this far'
!!$    stop

!!$  Now we actually do things

    t = 0.0_wp

    call setup_base_grid(base)

!!$    u   => base%u
!!$    u_p => base%u_p
!!$    x   => base%x

    allocate(u_exact(base%nvars, 1-ghostzones:base%n_points+ghostzones),&
             error  (base%nvars, 1-ghostzones:base%n_points+ghostzones),&
             STAT=ierr)

    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate arrays!'
      STOP
    end if

    call exact_solution(base%nvars, base%n_points, ghostzones, u_exact, &
         base%x, t)
    base%u   = u_exact
    call exact_solution(base%nvars, base%n_points, ghostzones, u_exact, &
         base%x, t - base%dt)
    base%u_p = u_exact
    
!!$    u   = u_exact
!!$    u_p = u_exact

!!$    Initialize child

    child_box%lbnd = (base%mybox%lbnd + base%mybox%ubnd + 1) / 2 - &
         base%n_points / 4
    child_box%ubnd = base%mybox%ubnd + base%n_points / 2
    child_box%npoints = base%mybox%npoints
    child_box%nghosts = base%mybox%nghosts
    child_box%bbox = 0
    parent => base

!!$    call setup_grid(child, child_box, &
!!$         base%nvars, &
!!$         base%dx / 2.0_wp, &
!!$         base%dt / 2.0_wp, &
!!$         parent)
    call setup_grid(child, child_box, parent)

    call prolongate_grid(child)

!!$  Having initialized, run until the end time 

    do while ( last_time - t > 0.000001_wp * base%dt)

      call swap_timelevels(base)

!!$      u   => base%u
!!$      u_p => base%u_p
      
      t = t + base%dt

!!$      Single evolution step

      call evolve_one_step(base)

!!$      Evolve child grid

      call swap_timelevels(child)
      call evolve_one_step(child)
      call swap_timelevels(child)
      call evolve_one_step(child)
      
      call restrict_grid(child)

!!$      Look at the AMR errors

!!$      call find_error_box(child, child_box)

    end do
      
    call exact_solution(base%nvars, base%n_points, ghostzones, u_exact, &
         base%x, t)

    error = abs(base%u - u_exact)

    write(*,*) 'Outputting for n = ', nconv

    call output('tvd_base', base%nvars, base%n_points, ghostzones, &
         base%x, base%u, u_exact, error)

    if (output_exact_and_errors) then

      call output_errors('tvd_base', base%nvars, base%n_points, ghostzones, &
           nconv, base%dx, error)

    end if
    
    call exact_solution(child%nvars, child%n_points, ghostzones, u_exact, &
         child%x, t)

    error = abs(child%u - u_exact)

    call output('tvd_child', child%nvars, child%n_points, ghostzones, &
         child%x, child%u, u_exact, error)

    if (output_exact_and_errors) then

      call output_errors('tvd_child', child%nvars, child%n_points, &
           ghostzones, &
           nconv, child%dx, error)

    end if

!!$    write(*,*) 'Deleting child grid'
    call delete_grid(child)
!!$    write(*,*) 'Deleting base grid'
    call delete_grid(base)

    deallocate(u_exact,&
               error)

  end do

end program tvd
