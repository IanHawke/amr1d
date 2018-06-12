!!$ This module defines a single grid.
!!$ It should contain the mesh refinement information, grid spacing etc.

module grid_mod

  use real_type_mod
  use grid_parameters_mod
  use system_parameters_mod
  use data_mod
  use box_mod

  type grid

     integer :: n_points
     integer :: ghostzones
     integer :: nvars

!!$iteration is the current iteration that this grid is on.
!!$iteration_step is the number of iterations that a single evolution step
!!$will increase iteration by, at this refinement level.

     integer :: iteration
     integer :: iteration_step

!!$x0 is the location of the left boundary of the domain covered by this grid

     real(kind=wp) :: x0

     real(kind=wp) :: dx, dt

     real(kind=wp), dimension(:,:), pointer :: u, u_p, local_error

     real(kind=wp), dimension(:), pointer :: x
              
     type(data) :: grid_data

     type(grid), pointer :: parent

     type(box) :: mybox

!!$     integer :: grid_id
                               
  end type grid

!!$  integer, save :: grid_id_max = 1

contains

  subroutine setup_base_grid(base)

    implicit none

    integer :: i, ierr

    type(grid), intent(out) :: base

    base%n_points       = base_n_points
    base%ghostzones     = ghostzones
    base%nvars          = nvariables
    base%iteration      = 0
    base%iteration_step = 2**(max_refinement_levels - 1)
    base%dt             = base_dt
    base%dx             = base_dx
    base%x0             = xmin

    base%mybox%lbnd    = 1
    base%mybox%ubnd    = base%n_points
    base%mybox%npoints = base%n_points
    base%mybox%nghosts = base%ghostzones
    base%mybox%bbox    = 1

!!$    base%grid_id       = grid_id_max
!!$    grid_id_max = grid_id_max + 1

    allocate(base%u  (base%nvars,1-base%ghostzones:base%n_points+base%ghostzones),&
             base%u_p(base%nvars,1-base%ghostzones:base%n_points+base%ghostzones),&
             base%local_error(base%nvars,1-base%ghostzones:base%n_points+base%ghostzones),&
             base%x  (1-base%ghostzones:base%n_points+base%ghostzones),&
             STAT=ierr)

    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate arrays!'
      STOP
    end if
    
    do i = 1 - base%ghostzones, base%n_points + base%ghostzones
      base%x(i) = base%x0 + (dble(i) - 0.5_wp) * base%dx
    end do

!!$Ensure that the base grid is always completely regridded for self-shadowing

    base%local_error = 10.0_wp * amr_error
    
    call setup_data(base%nvars, base%n_points, base%ghostzones, &
                    base%grid_data)

    nullify(base%parent)
  
  end subroutine setup_base_grid

  subroutine setup_grid(new_grid, new_box, parent)

    implicit none

    integer :: i, ierr

    type(grid), intent(out)   :: new_grid
    type(box), intent(in)     :: new_box
    type(grid), pointer       :: parent

    new_grid%n_points       = new_box%npoints
    new_grid%ghostzones     = new_box%nghosts
    new_grid%nvars          = parent%nvars
    new_grid%iteration      = parent%iteration
    if (parent%iteration_step == 1) then
      write(*,*) 'Tried to allocate a level finer than max_refinement_levels allows!'
      STOP
    end if
    new_grid%iteration_step = parent%iteration_step / 2
    new_grid%dt             = parent%dt / 2.0_wp
    new_grid%dx             = parent%dx / 2.0_wp
    new_grid%parent         => parent

    new_grid%mybox = new_box

!!$    new_grid%grid_id        = grid_id_max
!!$    grid_id_max = grid_id_max + 1
!!$    write(*,*) 'Set up grid id  ', new_grid%grid_id

    allocate(new_grid%u  (new_grid%nvars,1-new_grid%ghostzones:new_grid%n_points+new_grid%ghostzones),&
             new_grid%u_p(new_grid%nvars,1-new_grid%ghostzones:new_grid%n_points+new_grid%ghostzones),&
             new_grid%local_error(new_grid%nvars,1-new_grid%ghostzones:new_grid%n_points+new_grid%ghostzones),&
             new_grid%x  (1-new_grid%ghostzones:new_grid%n_points+new_grid%ghostzones),&
             STAT=ierr)

    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate arrays!'
      STOP
    end if

    new_grid%x0 = parent%x0 + (new_box%lbnd - 1) * parent%dx

    do i = 1 - new_grid%ghostzones, new_grid%n_points + new_grid%ghostzones
      new_grid%x(i) = new_grid%x0 + (dble(i) - 0.5_wp) * new_grid%dx
    end do

    call setup_data(new_grid%nvars, new_grid%n_points, new_grid%ghostzones, &
                    new_grid%grid_data)
  
  end subroutine setup_grid

  subroutine delete_grid(base)

    implicit none

    integer :: i, ierr
    
    type(grid), intent(out) :: base

    ierr = 0

    base%n_points   = -1
    base%ghostzones = -1
    base%dt         = -1.0_wp
    base%dx         = -1.0_wp

    base%parent     => NULL()
    
    if (associated(base%u, base%u_p)) then
      write(*,*) 'The two pointers are associated here!!!'
    end if
    
    if (.not.associated(base%u)) then
      write(*,*) 'base%u not associated'
    end if
    
    if (.not.associated(base%u_p)) then
      write(*,*) 'base%u_p not associated'
    end if
    
    deallocate(base%u          ,&
               base%u_p        ,&
               base%local_error,&
               base%x          ,&
               STAT=ierr)
    
    if (ierr .ne. 0) then
      write(*,*) 'Failed to deallocate arrays in delete_grid!'
      STOP
    end if
    
    call delete_data(base%grid_data)
    
  end subroutine delete_grid
  
  subroutine initialize_grid(base)

    implicit none
    
    type(grid) :: base

    real(kind=wp), allocatable, dimension(:,:) :: u_exact
    real(kind=wp) :: t = 0.0_wp
    integer :: ierr

    allocate(u_exact(base%nvars, &
                     1 - base%ghostzones:base%n_points + base%ghostzones),&
             STAT = ierr)

    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate arrays (in initialize_base_grid)!'
      STOP
    end if

    call exact_solution(base%nvars, base%n_points, base%ghostzones, u_exact, &
         base%x, t)
    base%u   = u_exact
    call exact_solution(base%nvars, base%n_points, base%ghostzones, u_exact, &
         base%x, t - base%dt)
    base%u_p = u_exact

    deallocate(u_exact)

  end subroutine initialize_grid
  
  subroutine swap_timelevels(base)
    
    implicit none
    
    type(grid) :: base
    
    real(kind=wp), dimension(:,:), pointer :: u_tmp
    
    u_tmp    => base%u
    base%u   => base%u_p
    base%u_p => u_tmp
    
  end subroutine swap_timelevels

!!$Prolongation and restriction operators

  subroutine prolongate_grid(child)

    implicit none

    type(grid), intent(inout) :: child

    real(kind=wp) :: slope_l, slope_r, slope

    integer :: i, k, parent_i, istart, iend

    if (mod(child%n_points, 2) .ne. 0) then
      write(*,*) 'Child grids must be an even number of points!'
      write(*,*) 'That is, every parent cell must either be totally refined or not'
      STOP
    end if
    
!!$    Piecewise constant prolongation

!!$    do i = 1 - child%ghostzones, child%n_points + child%ghostzones, 2
!!$
!!$      child%u(:, i    ) = child%parent%u(:, child%mybox%lbnd + (i - 1) / 2)
!!$      child%u(:, i + 1) = child%parent%u(:, child%mybox%lbnd + (i - 1) / 2)
!!$
!!$    end do

!!$    Linear prolongation using minmod

    if (mod(child%ghostzones, 2) == 0) then

       istart = 1 - child%ghostzones
       iend   = child%n_points + child%ghostzones

    else

       istart = 1 - child%ghostzones + 1
       iend   = child%n_points + child%ghostzones - 1

       do k = 1, child%nvars

!!$Left most ghostzone

          parent_i = child%mybox%lbnd - (child%ghostzones + 1) / 2

          slope_l = &
               child%parent%u(k, parent_i    ) - child%parent%u(k, parent_i - 1)
          slope_r = &
               child%parent%u(k, parent_i + 1) - child%parent%u(k, parent_i    )
          
          if (slope_l * slope_r < 0.0_wp) then
             slope = 0.0_wp
          else if(abs(slope_l) < abs(slope_r)) then
             slope = slope_l
          else
             slope = slope_r
          end if
          
          child%u(k, 1 - child%ghostzones)   = child%parent%u(k, parent_i) + &
               0.25_wp * slope

!!$Right most ghostzone

          parent_i = child%mybox%lbnd +child%n_points + child%ghostzones / 2

          slope_l = &
               child%parent%u(k, parent_i    ) - child%parent%u(k, parent_i - 1)
          slope_r = &
               child%parent%u(k, parent_i + 1) - child%parent%u(k, parent_i    )
          
          if (slope_l * slope_r < 0.0_wp) then
             slope = 0.0_wp
          else if(abs(slope_l) < abs(slope_r)) then
             slope = slope_l
          else
             slope = slope_r
          end if
          
          child%u(k, child%n_points + child%ghostzones)   = &
               child%parent%u(k, parent_i) - 0.25_wp * slope
       end do
       
    end if

!!$Deal with all points (except for final odd numbered ghostzones, done above
    
    do i = istart, iend, 2
      do k = 1, child%nvars

        parent_i = child%mybox%lbnd + (i - 1) / 2

        slope_l = &
             child%parent%u(k, parent_i    ) - child%parent%u(k, parent_i - 1)
        slope_r = &
             child%parent%u(k, parent_i + 1) - child%parent%u(k, parent_i    )

        if (slope_l * slope_r < 0.0_wp) then
          slope = 0.0_wp
        else if(abs(slope_l) < abs(slope_r)) then
          slope = slope_l
        else
          slope = slope_r
        end if
        
        child%u(k, i)       = child%parent%u(k, parent_i) - &
             0.25_wp * slope
        child%u(k, i + 1)   = child%parent%u(k, parent_i) + &
             0.25_wp * slope

      end do
    end do
    
  end subroutine prolongate_grid

!!$ Prolong just the boundaries. Linear in time only.
!!$ This routine takes an extra argument. The time of the child grid
!!$ is not completely specified by the grid iteration if it is in the
!!$ middle of the RK step. Instead, it may be offset by
!!$ tplus \in [0,1], so
!!$ child time = child_iteration_time + child_dt * tplus, in effect.

  subroutine prolongate_boundaries(child, child_data, tplus)

    implicit none

    type(grid), intent(in)       :: child
    real(kind=wp), &
         dimension(child%nvars, &
                   1 - child%ghostzones : child%n_points + child%ghostzones), &
         intent(inout)           :: child_data
    real(kind=wp), intent(inout) :: tplus

    real(kind=wp), parameter :: grid_tiny = 0.0000000000001_wp

    integer :: i, k, parent_i, istart, iend
    real(kind=wp) :: u_weight, u_p_weight
    real(kind=wp) :: slope_l, slope_r, slope
    real(kind=wp) :: child_u_i, child_u_p_i, child_u_ip1, child_u_p_ip1

    logical :: do_left, do_right

!!$    Do some checking

    if (.not.associated(child%parent)) then
      write(*,*) 'This child does not have an associated parent!'
      STOP
    end if
    
    if (child%parent%dx < 0.0_wp) then
      write(*,*) 'This child has a deleted parent!!'
      STOP
    end if
    
    if (child%iteration > child%parent%iteration) then
      write(*,*) 'Child is ahead of parent in boundary prolongation!'
      STOP
    end if

    if (child%iteration < &
        child%parent%iteration - child%parent%iteration_step) then
      write(*,*) &
           'Child is not between parent time steps in boundary prolongation!'
      STOP
    end if
    
    if ((tplus < -grid_tiny) .or. (tplus > 1.0 + grid_tiny)) then
      write(*,*) 'tplus out of range in boundary prolongation!'
      STOP
    end if
    
    if (tplus < 0.0_wp) then
      tplus = 0.0_wp
    end if
    if (tplus > 1.0_wp) then
      tplus = 1.0_wp
    end if

!!$    Check if a boundary is an outer boundary

    do_left = .true.
    if (child%mybox%bbox(1) == 1) then
      do_left = .false.
    end if
    do_right = .true.
    if (child%mybox%bbox(2) == 1) then
      do_right = .false.
    end if
    
!!$    Find the weights for the interpolation in time

    if (child%iteration == child%parent%iteration) then
      
      if (abs(tplus) > grid_tiny) then
        write(*,*) 'Parent and child aligned in time but tplus not zero!'
        STOP
      end if
      
      u_weight   = 1.0_wp
      u_p_weight = 0.0_wp

    else if (child%iteration == &
         child%parent%iteration - child%iteration_step) then

      u_weight   = 0.5_wp * (1.0_wp + tplus)
      u_p_weight = 0.5_wp * (1.0_wp - tplus)

    else if (child%iteration == &
         child%parent%iteration - child%parent%iteration_step) then

      u_weight   = 0.5_wp * tplus
      u_p_weight = 0.5_wp * (2.0_wp - tplus)

    else

      write(*,*) 'Self checks failed in boundary prolongation!'
      STOP

    end if

    if (mod(child%ghostzones, 2) == 0) then

       istart = 1
       iend   = child%ghostzones

    else

       istart = 1
       iend   = child%ghostzones - 1

       do k = 1, child%nvars

!!$Left most ghostzone

         if (do_left) then

           parent_i = child%mybox%lbnd - (child%ghostzones + 1) / 2

           slope_l = &
                child%parent%u(k, parent_i    ) - child%parent%u(k, parent_i - 1)
           slope_r = &
                child%parent%u(k, parent_i + 1) - child%parent%u(k, parent_i    )
          
           if (slope_l * slope_r < 0.0_wp) then
             slope = 0.0_wp
           else if(abs(slope_l) < abs(slope_r)) then
             slope = slope_l
           else
             slope = slope_r
           end if
           
           child_u_ip1 = child%parent%u(k, parent_i) + &
                0.25_wp * slope
           
           slope_l = &
                child%parent%u_p(k, parent_i    ) - child%parent%u_p(k, parent_i - 1)
           slope_r = &
                child%parent%u_p(k, parent_i + 1) - child%parent%u_p(k, parent_i    )
           
           if (slope_l * slope_r < 0.0_wp) then
             slope = 0.0_wp
           else if(abs(slope_l) < abs(slope_r)) then
             slope = slope_l
           else
             slope = slope_r
           end if
           
           child_u_p_ip1 = child%parent%u_p(k, parent_i) + &
                0.25_wp * slope
           
           child_data(k, 1 - child%ghostzones)   = u_weight   * child_u_ip1 + &
                                                   u_p_weight * child_u_p_ip1

         end if
         
!!$Right most ghostzone

         if (do_right) then

           parent_i = child%mybox%lbnd + &
                (child%n_points + child%ghostzones) / 2

           slope_l = &
                child%parent%u(k, parent_i    ) - child%parent%u(k, parent_i - 1)
           slope_r = &
                child%parent%u(k, parent_i + 1) - child%parent%u(k, parent_i    )
          
           if (slope_l * slope_r < 0.0_wp) then
             slope = 0.0_wp
           else if(abs(slope_l) < abs(slope_r)) then
             slope = slope_l
           else
             slope = slope_r
           end if
           
           child_u_i = &
                child%parent%u(k, parent_i) - 0.25_wp * slope

           slope_l = &
                child%parent%u_p(k, parent_i    ) - child%parent%u_p(k, parent_i - 1)
           slope_r = &
                child%parent%u_p(k, parent_i + 1) - child%parent%u_p(k, parent_i    )
          
           if (slope_l * slope_r < 0.0_wp) then
             slope = 0.0_wp
           else if(abs(slope_l) < abs(slope_r)) then
             slope = slope_l
           else
             slope = slope_r
           end if
           
           child_u_p_i = &
                child%parent%u_p(k, parent_i) - 0.25_wp * slope

           child_data(k, child%n_points + child%ghostzones)   = &
                u_weight   * child_u_i + &
                u_p_weight * child_u_p_i

         end if
         
       end do
       
    end if

!!$    Prolongation at the left boundary
    
    do i = istart, iend, 2
      do k = 1, child%nvars

!!$         Left ghost zones

        if (do_left) then

          parent_i = child%mybox%lbnd - 1 - (i - 1) / 2

          slope_l = &
               child%parent%u(k, parent_i    ) - child%parent%u(k, parent_i - 1)
          slope_r = &
               child%parent%u(k, parent_i + 1) - child%parent%u(k, parent_i    )
          
          if (slope_l * slope_r < 0.0_wp) then
            slope = 0.0_wp
          else if(abs(slope_l) < abs(slope_r)) then
            slope = slope_l
          else
            slope = slope_r
          end if
          
          child_u_i   = child%parent%u(k, parent_i) - &
               0.25_wp * slope
          child_u_ip1 = child%parent%u(k, parent_i) + &
               0.25_wp * slope
          
          slope_l = &
               child%parent%u_p(k, parent_i    ) - &
               child%parent%u_p(k, parent_i - 1)
          slope_r = &
               child%parent%u_p(k, parent_i + 1) - &
               child%parent%u_p(k, parent_i    )
          
          if (slope_l * slope_r < 0.0_wp) then
            slope = 0.0_wp
          else if(abs(slope_l) < abs(slope_r)) then
            slope = slope_l
          else
            slope = slope_r
          end if
          
          child_u_p_i   = child%parent%u_p(k, parent_i) - &
               0.25_wp * slope
          child_u_p_ip1 = child%parent%u_p(k, parent_i) + &
               0.25_wp * slope
          
          child_data(k, 0 - i) = u_weight   * child_u_i + &
                                 u_p_weight * child_u_p_i
          child_data(k, 1 - i) = u_weight   * child_u_ip1 + &
                                 u_p_weight * child_u_p_ip1

        end if
        
!!$         Right ghost zones

        if (do_right) then

          parent_i = child%mybox%lbnd + (child%n_points + i) / 2

          slope_l = &
               child%parent%u(k, parent_i    ) - child%parent%u(k, parent_i - 1)
          slope_r = &
               child%parent%u(k, parent_i + 1) - child%parent%u(k, parent_i    )
          
          if (slope_l * slope_r < 0.0_wp) then
            slope = 0.0_wp
          else if(abs(slope_l) < abs(slope_r)) then
            slope = slope_l
          else
            slope = slope_r
          end if
          
          child_u_i   = child%parent%u(k, parent_i) - &
               0.25_wp * slope
          child_u_ip1 = child%parent%u(k, parent_i) + &
               0.25_wp * slope
          
          slope_l = &
               child%parent%u_p(k, parent_i    ) - &
               child%parent%u_p(k, parent_i - 1)
          slope_r = &
               child%parent%u_p(k, parent_i + 1) - &
               child%parent%u_p(k, parent_i    )
          
          if (slope_l * slope_r < 0.0_wp) then
            slope = 0.0_wp
          else if(abs(slope_l) < abs(slope_r)) then
            slope = slope_l
          else
            slope = slope_r
          end if
          
          child_u_p_i   = child%parent%u_p(k, parent_i) - &
               0.25_wp * slope
          child_u_p_ip1 = child%parent%u_p(k, parent_i) + &
               0.25_wp * slope

          child_data(k, child%n_points + i)     = u_weight   * child_u_i + &
                                                  u_p_weight * child_u_p_i
          child_data(k, child%n_points + i + 1) = u_weight   * child_u_ip1 + &
                                                  u_p_weight * child_u_p_ip1
        
        end if
        
      end do
    end do
    
  end subroutine prolongate_boundaries

!!$Restriction operator

  subroutine restrict_grid(child)

    implicit none

    type(grid), intent(inout) :: child

    integer :: i, k, parent_i, istart, iend

    if (mod(child%n_points, 2) .ne. 0) then
       write(*,*) 'The child grid should always have an even number of points!'
       STOP
    end if
    
    do i = 1, child%n_points, 2

      parent_i = child%mybox%lbnd + (i - 1) / 2

      do k = 1, child%nvars

        child%local_error(k, i) = child%parent%u(k, parent_i)

        child%parent%u(k, parent_i) = 0.5_wp * ( &
             child%u(k, i    ) + &
             child%u(k, i + 1)                 )

        child%local_error(k, i) = &
             abs( child%parent%u   (k, parent_i) - &
                  child%local_error(k, i       ) )

        child%local_error(k, i + 1) = child%local_error(k, i)
        
      end do
      
    end do
    
  end subroutine restrict_grid

!!$From the local errors work out where the grid should be refined

  subroutine find_error_boxes(base, child_box_list)

    use box_list_mod

    implicit none

    type(grid), intent(in)      :: base
    type(list_type), intent(out) :: child_box_list

    type(box), pointer :: child_box
  
    logical, dimension(:), allocatable :: in_error, error_tmp

    integer :: ierr
    integer :: i, k, istart, iend
    real(kind=wp) :: point_error

    allocate(in_error(1 - base%ghostzones: base%n_points + base%ghostzones), &
            error_tmp(1 - base%ghostzones: base%n_points + base%ghostzones), &
            STAT = ierr)

    if (ierr .ne. 0) then
       write(*,*) 'Error allocating arrays in find_error_box', ierr
       STOP
    end if

    in_error = .false.
    error_tmp = .false.

!!$For the moment this is just absolute error. Needs generalizing

    do i = 1, base%n_points

       point_error = 0.0_wp
       do k = 1, base%nvars
!!$         if (((euler).and.(k==3)).or.(.not.euler)) then
           point_error = point_error + base%local_error(k, i)
!!$         end if
       end do
       
       if (point_error > amr_error) then
          error_tmp(i) = .true.
       end if

    end do
    
!!$Buffering

    in_error = error_tmp
    
    do i = 1, base%n_points

       istart = max(1, i - amr_buffer_width)
       iend   = min(base%n_points, i + amr_buffer_width)

       if (any(error_tmp(istart:iend))) then
          in_error(i) = .true.
       end if
       
    end do
    
!!$Boundary check - needs extending for multiple models

!!$    if (periodic) then
!!$      if (any(in_error(1:amr_buffer_width)).or.&
!!$          any(in_error(base%n_points - amr_buffer_width + 1: base%n_points))) &
!!$        then
!!$          in_error(1:amr_buffer_width) &
!!$               = .true.
!!$          in_error(base%n_points - amr_buffer_width + 1: base%n_points) & 
!!$               = .true.
!!$      end if
!!$    end if

!!$    do i = 1, base%n_points
!!$
!!$       write(15, *) base%x(i), base%local_error(1, i), error_tmp(i), in_error(i)
!!$       
!!$    end do
!!$
!!$    write(15, *)
   
!!$Set the boxes 

    istart = 0
    iend   = -1

    do i = 1, base%n_points
!!$       if we have not found the start of the box, search
      if (istart < 1) then
        if (in_error(i)) then
          istart = i
        end if
      else
!!$       we have found the start - need to find the finish
        if (.not.in_error(i)) then
          iend = i - 1
!!$             Now we need to compute the box etc.
          allocate(child_box)
          child_box%lbnd = istart
          child_box%ubnd = iend
          child_box%npoints = 2*(iend - istart + 1)
          child_box%nghosts = base%ghostzones
          if (istart == 1) then
            child_box%bbox(1) = base%mybox%bbox(1)
          else
            child_box%bbox(1) = 0
          end if
!!$             If we have actually found the end then this cannot
!!$             be true, but we leave it anyway.
          if (iend == base%n_points) then
            child_box%bbox(2) = base%mybox%bbox(2)
          else
            child_box%bbox(2) = 0
          end if
          call add_box_to_list(child_box_list, child_box)
          child_box => NULL()
          istart = 0
          iend = -1
        end if
      end if
    end do
     
!!$    Check that there is not a box to the right boundary
    if (istart > 0) then
      iend = base%n_points
!!$             Now we need to compute the box etc.
      allocate(child_box)
      child_box%lbnd = istart
      child_box%ubnd = iend
      child_box%npoints = 2*(iend - istart + 1)
      child_box%nghosts = base%ghostzones
      if (istart == 1) then
        child_box%bbox(1) = base%mybox%bbox(1)
      else
        child_box%bbox(1) = 0
      end if
!!$      This should now always be true but etc.
      if (iend == base%n_points) then
        child_box%bbox(2) = base%mybox%bbox(2)
      else
        child_box%bbox(2) = 0
      end if
      call add_box_to_list(child_box_list, child_box)
      child_box => NULL()
      istart = 0
      iend = -1
    end if

    deallocate(in_error, error_tmp)
  
  end subroutine find_error_boxes

  subroutine output_grid(base, depth, ngrid, t)

    implicit none

    type(grid) :: base
    integer, intent(in) :: depth, ngrid
    real(kind=wp), intent(in) :: t

    real(kind=wp), allocatable, dimension(:,:) :: u_exact
    real(kind=wp), allocatable, dimension(:,:) :: error
    integer :: ierr

    character(len=200) :: filestring, &
         filestring_depth, filestring_depth_trim, &
         filestring_ngrid, filestring_ngrid_trim
    integer :: filestring_depth_len, filestring_ngrid_len

    allocate(u_exact(base%nvars, 1-base%ghostzones:base%n_points+base%ghostzones),&
             error  (base%nvars, 1-base%ghostzones:base%n_points+base%ghostzones),&
             STAT=ierr)

    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate arrays!'
      STOP
    end if

    call exact_solution(base%nvars, base%n_points, base%ghostzones, u_exact, &
         base%x, t)

    error = abs(base%u - u_exact)

    write(filestring_depth, *) depth
    filestring_depth_len = int(log10(dble(depth)))+1
    filestring_depth_trim = trim( &
       filestring_depth(len_trim(filestring_depth) + 1 - filestring_depth_len:&
       len_trim(filestring_depth)))
    write(filestring_ngrid, *) ngrid
    filestring_ngrid_len = int(log10(dble(ngrid)))+1
    filestring_ngrid_trim = trim( &
       filestring_ngrid(len_trim(filestring_ngrid) + 1 - filestring_ngrid_len:&
       len_trim(filestring_ngrid)))

    filestring = trim('amr_l'//filestring_depth_trim(1:filestring_depth_len)//&
                      '_g'//filestring_ngrid_trim(1:filestring_ngrid_len)//'_')

    call output(trim(filestring), base%nvars, base%n_points, base%ghostzones, &
         base%x, base%u, u_exact, error)

    deallocate(u_exact, error)

  end subroutine output_grid

!!$  This checks if the passed in "parent" would be a valid parent

  logical function grid_child_of_parent(child, parent)

    implicit none

    type(grid) :: child, parent

    grid_child_of_parent = .false.

    if ((child%x0 - 0.5_wp * child%dx > parent%x0 - 0.5_wp * parent%dx).and.&
        (child%x0  + child%n_points  * child%dx  + 0.5_wp * child%dx < &
         parent%x0 + parent%n_points * parent%dx + 0.5_wp * parent%dx)) then
      grid_child_of_parent = .true.
    end if
    
  end function grid_child_of_parent

!!$This child grid has the correct parent grid, but the box entry may
!!$be wrong after reconnection.

  subroutine correct_box(child)

    implicit none

    type(grid) :: child

    integer :: offset

    offset = nint((child%x0 - child%parent%x0) / child%parent%dx)
    child%mybox%lbnd = offset + 1
    child%mybox%ubnd = child%mybox%lbnd + child%n_points / 2 - 1
  
  end subroutine correct_box

  subroutine set_initial_error_on_grid(base)

    implicit none

    type(grid), pointer :: base

    integer :: i

    do i = 2, base%n_points

      base%local_error(:, i) = abs(base%u(:, i) - base%u(:, i - 1))

    end do
    
  end subroutine set_initial_error_on_grid

  subroutine apply_physical_boundaries_grid(base)

    implicit none

    type(grid) :: base

    if ((base%mybox%bbox(1) .ne. 0) .and.&
        (base%mybox%bbox(2) .ne. 0)) then
      call boundary_conditions        (base%nvars, &
                                       base%n_points, &
                                       base%ghostzones, &
                                       base%u)
    else
      if (base%mybox%bbox(1) .ne. 0) then
        call boundary_conditions_left (base%nvars, &
                                       base%n_points, &
                                       base%ghostzones, &
                                       base%u)
      end if
      if (base%mybox%bbox(2) .ne. 0) then
        call boundary_conditions_right(base%nvars, &
                                       base%n_points, &
                                       base%ghostzones, &
                                       base%u)
      end if
    end if
    

  end subroutine apply_physical_boundaries_grid

!!$  This routine looks for the intersection between two grids at the same
!!$  depth. It returns the start and end locations inside the "to" grid of
!!$  the intersection, and the offset.

  logical function grid_intersection(grid_to, grid_from, istart, iend, offset)

    implicit none

    type(grid), pointer :: grid_to, grid_from
    integer, intent(out) :: istart, iend, offset

    real(kind=wp) :: to_start, to_end, from_start, from_end, dx
    real(kind=wp) :: intersection_start, intersection_end

    istart = -10
    iend   = -10
    offset = 0

    grid_intersection = .false.

    dx = grid_to%dx
    to_start   = grid_to%x0
    to_end     = grid_to%x0   + grid_to%n_points   * dx
    from_start = grid_from%x0
    from_end   = grid_from%x0 + grid_from%n_points * dx

    if (to_start + 0.5_wp * dx < from_end) then
      if (to_end > from_start + 0.5_wp * dx) then
        grid_intersection = .true.
        intersection_start = max(to_start, from_start)
        intersection_end   = min(to_end  , from_end  )
        istart = nint((intersection_start - to_start) / dx) + 1
        offset = nint((to_start - from_start) / dx)
        iend = nint((intersection_end - to_start) / dx)
      end if
    end if
    
  end function grid_intersection

  subroutine set_dt_from_parent_grid(child_grid)

    implicit none

    type(grid) :: child_grid

    child_grid%dt = child_grid%parent%dt / 2.0_wp

  end subroutine set_dt_from_parent_grid

  subroutine find_max_char_grid(base, max_char)
  
    implicit none

    type(grid) :: base
    real(kind=wp), intent(inout) :: max_char

    integer :: nv, np, ngz
    real(kind=wp) :: char_max

    nv  = base%nvars
    np  = base%n_points
    ngz = base%ghostzones

    call maxval_char(nv, np, ngz, base%u, char_max)
      
    max_char = max(max_char, char_max)

  end subroutine find_max_char_grid
  
end module grid_mod
