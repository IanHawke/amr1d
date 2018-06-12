!!$ This module defines data (utility GFs that have only one timelevel).
!!$ Each grid should have one.

module data_mod

  use real_type_mod
  use method_parameters_mod

  type data


     real(kind=wp), dimension(:,:), pointer :: u_int1, u_int2, u_int3, &
                                               rhs, &
                                               flux, &
                                               f_plus, f_minus, &
                                               f_plus_right, f_minus_left, &
                                               f_char_plus_right, f_char_minus_left, &
                                               ul, ur, &
                                               slope, &
                                               slope_l, slope_r

     real(kind=wp), dimension(:,:,:), pointer :: u_char, f_char, &
                                                 f_char_minus, f_char_plus


     real(kind=wp), dimension(:,:,:), pointer :: leftev, rightev

  end type data

contains

  subroutine setup_data(nv, np, ngz, base)

    implicit none
    
    integer :: ierr
    
    integer, intent(in) :: nv, np, ngz
    type(data), intent(out) :: base
    
    
    if (RK2) then
      
      allocate(base%u_int1(nv,1-ngz:np+ngz),&
               STAT=ierr)
      
      if (ierr .ne. 0) then
        write(*,*) 'Failed to allocate arrays!'
        STOP
      end if
      
      nullify(base%u_int2)
      nullify(base%u_int3)
      
    else if (RK4) then
      
      allocate(base%u_int1(nv,1-ngz:np+ngz),&
               base%u_int2(nv,1-ngz:np+ngz),&
               base%u_int3(nv,1-ngz:np+ngz),&
               STAT=ierr)

      if (ierr .ne. 0) then
        write(*,*) 'Failed to allocate arrays!'
        STOP
      end if

    else

      write(*,*) 'How are you integrating??'
      STOP

    end if
      
    allocate(base%rhs (nv,1-ngz:np+ngz),&
             base%flux(nv,1-ngz:np+ngz),&
             base%ul  (nv,1-ngz:np+ngz),&
             base%ur  (nv,1-ngz:np+ngz),&
             STAT=ierr)
    
    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate arrays!'
      STOP
    end if

    if (TVD_differencing) then

      allocate(base%slope  (nv,1-ngz:np+ngz),&
               base%slope_l(nv,1-ngz:np+ngz),&
               base%slope_r(nv,1-ngz:np+ngz),&
               STAT=ierr)

      if (ierr .ne. 0) then
        write(*,*) 'Failed to allocate arrays!'
        STOP
      end if

    else
        
      nullify(base%slope)
      nullify(base%slope_l)
      nullify(base%slope_r)

    end if

    if (Finite_Difference) then

      allocate(base%f_minus     (nv,1-ngz:np+ngz),&
               base%f_plus      (nv,1-ngz:np+ngz),&
               base%f_minus_left(nv,1-ngz:np+ngz),&
               base%f_plus_right(nv,1-ngz:np+ngz),&
               base%u_char(nv,-2:3,1-ngz:np+ngz),&
               base%f_char(nv,-2:3,1-ngz:np+ngz),&
               base%f_char_minus(nv,-2:3,1-ngz:np+ngz),&
               base%f_char_plus (nv,-2:3,1-ngz:np+ngz),&
               base%f_char_minus_left(nv,1-ngz:np+ngz),&
               base%f_char_plus_right(nv,1-ngz:np+ngz),&
               base%leftev (nv,nv,1-ngz:np+ngz),&
               base%rightev(nv,nv,1-ngz:np+ngz),&
               STAT=ierr)

      if (ierr .ne. 0) then
        write(*,*) 'Failed to allocate arrays!'
        STOP
      end if
      
    else
        
      nullify(base%f_minus)
      nullify(base%f_plus )
      nullify(base%f_minus_left)
      nullify(base%f_plus_right)
      nullify(base%u_char)
      nullify(base%f_char)
      nullify(base%f_char_minus)
      nullify(base%f_char_plus)
      nullify(base%f_char_minus_left)
      nullify(base%f_char_plus_right)
      nullify(base%leftev)
      nullify(base%rightev)

    end if
      
  end subroutine setup_data
    
  subroutine delete_data(base)

    implicit none

    integer :: ierr

    type(data), intent(out) :: base
      
    if (RK2) then

      deallocate(base%u_int1)

    else if (RK4) then
      
      deallocate(base%u_int1, &
                 base%u_int2, &
                 base%u_int3)

    else

      write(*,*) 'How are you integrating??'
      STOP

    end if

    deallocate(base%rhs ,&
               base%flux,&
               base%ul  ,&
               base%ur  ,&
               STAT=ierr)

    if (ierr .ne. 0) then
      write(*,*) 'Failed to deallocate arrays!'
      STOP
    end if

    if (TVD_differencing) then

      deallocate(base%slope  ,&
                 base%slope_l,&
                 base%slope_r,&
                 STAT=ierr)

      if (ierr .ne. 0) then
        write(*,*) 'Failed to deallocate arrays!'
        STOP
      end if

    end if

    if (Finite_Difference) then

      deallocate(base%f_plus  ,&
                 base%f_minus,&
                 base%f_plus_right,&
                 base%f_minus_left,&
                 base%u_char,&
                 base%f_char,&
                 base%f_char_minus,&
                 base%f_char_plus,&
                 base%f_char_minus_left,&
                 base%f_char_plus_right,&
                 base%leftev,&
                 base%rightev,&
                 STAT=ierr)

      if (ierr .ne. 0) then
        write(*,*) 'Failed to deallocate arrays!'
        STOP
      end if

    end if

  end subroutine delete_data
    
end module data_mod
