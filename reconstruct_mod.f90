!!$Reconstruct the function to the cell boundaries

module reconstruct_mod

contains

  subroutine reconstruct(nv, np, ngz, u, grid_data)

    use real_type_mod
    use method_parameters_mod
    use data_mod
    
    implicit none
    
    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: u
    type(data) :: grid_data
    
    integer :: ierr
    
    if (TVD_Reconstruct) then
      
      call find_slopes(nv, np, ngz, u, grid_data%ul, grid_data%ur, &
                       grid_data%slope, &
                       grid_data%slope_l, &
                       grid_data%slope_r)

    else if (WENO_Reconstruct) then

      call weno_reconstruction(nv, np, ngz, u, grid_data%ul, grid_data%ur)
    
    else
      
      write(*,*) "Reconstruction method not recognized"
      STOP
      
    end if

    if (reconstruct_primitive_vars) then

      call primitive2conservative(nv, np, ngz, grid_data%ul)
      call primitive2conservative(nv, np, ngz, grid_data%ur)
      
    else
      
      call conservative2primitive(nv, np, ngz, grid_data%ul)
      call conservative2primitive(nv, np, ngz, grid_data%ur)
      
    end if
    
  end subroutine reconstruct

  subroutine reconstruct_split_flux(nv, np, ngz, grid_data)

    use real_type_mod
    use method_parameters_mod
    use data_mod
  
    implicit none
    
    integer, intent(in) :: nv, np, ngz
    type(data) :: grid_data
    
    real(kind=wp), dimension(:, :), allocatable :: tmp

    integer :: ierr, var, i

    allocate(tmp(nv, 1 - ngz:np + ngz), &
             STAT=ierr)
    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate array in reconstruct_split_flux'
      STOP
    end if

    if (TVD_Reconstruct) then

      call find_slopes(nv, np, ngz, grid_data%f_minus, &
                       grid_data%f_minus_left, tmp, &
                       grid_data%slope, grid_data%slope_l, grid_data%slope_r)

      call find_slopes(nv, np, ngz, grid_data%f_plus, &
                       tmp, grid_data%f_plus_right, &
                       grid_data%slope, grid_data%slope_l, grid_data%slope_r)

    else if (WENO_Reconstruct) then

!!$      This is extremely dangerous - really we should be reconstructing
!!$      the characteristic fluxes
!!$      We are now doing this at the higher level following Shu

      STOP

!!$      grid_data%fr_char = 0.0_wp
!!$      grid_data%fl_char = 0.0_wp
!!$
!!$      do i = 0, np
!!$        do var = 1, nv
!!$          call WENO3_Right(-1, 3, &
!!$               grid_data%f_char_plus (var, -2:2, i), grid_data%fr_char(var, i-2:i+2))
!!$          call WENO3_Left (-1, 3, &
!!$               grid_data%f_char_minus(var, -1:3, i), grid_data%fl_char(var, i-1:i+3))
!!$        end do
!!$      end do

    else

      write(*,*) 'Do not recognize reconstruction method in reconstruct_split_flux'
      STOP

    end if
    
    deallocate(tmp)

  end subroutine reconstruct_split_flux
    
!!$This subroutine is added to dereference the pointers - it is purely
!!$a speed issue

  subroutine find_slopes(nv, np, ngz, u, ul, ur, slope, slope_l, slope_r)

    use real_type_mod
    use method_parameters_mod
    use data_mod
    
    implicit none
    
    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: u
    real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(out) :: ul, ur, &
                                                               slope, &
                                                               slope_l, &
                                                               slope_r

    integer :: ierr

!!$  Compute the slopes

    slope_l = 0.0_wp
    slope_r = 0.0_wp

    slope_l(:,0:np+1) = u(:,0:np+1) - u(:,-1:np  )
    slope_r(:,0:np+1) = u(:,1:np+2) - u(:, 0:np+1)
  
!!$  Limiters

    if (limiter_minmod) then

!!$    Minmod
  
      where (slope_l * slope_r < 0.0_wp)
        slope = 0.0_wp
      elsewhere (abs(slope_l) < abs(slope_r))
        slope = slope_l
      elsewhere
        slope = slope_r
      end where
  
    else if (limiter_mc) then

!!$    Van Leer MC limiter

      where (slope_l * slope_r < 0.0_wp)

        slope = 0.0_wp

      elsewhere

        slope = sign( min( 2.0_wp * abs(slope_l), &
                           2.0_wp * abs(slope_r), &
                           0.5_wp * (abs(slope_l) + abs(slope_r)) ), &
                         slope_l + slope_r ) 

      end where

    else
      
      write(*,*) 'Should never get to here!!!'
      
    end if
  
!!$  Reconstructed values

    ul = u - 0.5_wp * slope
    ur = u + 0.5_wp * slope
      
  end subroutine find_slopes

!!$WENO reconstruction - componentwise

  subroutine weno_reconstruction(nv, np, ngz, u, ul, ur)

    use real_type_mod
    use method_parameters_mod
    use data_mod
    
    implicit none
    
    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: u
    real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(out) :: ul, ur

    integer :: ierr, var

!!$Note the order seems odd - again

    do var = 1, nv
        
      call WENO3_Left (np, ngz, u(var, :), ur(var, :))
      call WENO3_Right(np, ngz, u(var, :), ul(var, :))

    end do
    
  end subroutine weno_reconstruction
  
  subroutine WENO3_Left(np, ngz, u, uminus)
    
    use real_type_mod
    
    implicit none
    
    integer :: i, np, ngz
    real(kind=wp), dimension(1 - ngz : np + ngz) :: u, uminus
    real(kind=wp) :: one, two, three, four, five, six, seven, &
         ten, eleven, twelve, thirteen, &
         ThirteenByTwelve, Quarter
    parameter (one = 1)
    parameter (two = 2)
    parameter (three = 3)
    parameter (four = 4)
    parameter (five = 5)
    parameter (six = 6)
    parameter (seven = 7)
    parameter (ten = 10)
    parameter (eleven = 11)
    parameter (twelve = 12)
    parameter (thirteen = 13)
    parameter (ThirteenByTwelve = thirteen / twelve)
    parameter (Quarter = one / four)
    real(kind=wp) :: d0, d1, d2
    parameter (d0 = three / ten)
    parameter (d1 = six / ten)
    parameter (d2 = one / ten)
    real(kind=wp) :: c00,c01,c02,c10,c11,c12,c20,c21,c22
    parameter (c00 = two / six)
    parameter (c01 = five / six)
    parameter (c02 = -one / six)
    parameter (c10 = -one / six)
    parameter (c11 = five / six)
    parameter (c12 = two / six)
    parameter (c20 = two / six)
    parameter (c21 = -seven / six)
    parameter (c22 = eleven / six)
    
    real(kind=wp) :: beta0, beta1, beta2
    real(kind=wp) :: epsilon
    real(kind=wp) :: alpha0, alpha1, alpha2, alphasum
    real(kind=wp) :: w0, w1, w2
    real(kind=wp) :: u0plushalf, u1plushalf, u2plushalf
    
    epsilon = 1.d-6
    
    do i = 0, np+1
      
      beta0 = ThirteenByTwelve * (u(i  ) - two * u(i+1) + u(i+2))**2 + &
           Quarter * (three * u(i  ) - four * u(i+1) +         u(i+2))**2
      beta1 = ThirteenByTwelve * (u(i-1) - two * u(i  ) + u(i+1))**2 + &
           Quarter * (        u(i-1)                 -         u(i+1))**2
      beta2 = ThirteenByTwelve * (u(i-2) - two * u(i-1) + u(i  ))**2 + &
           Quarter * (        u(i-2) - four * u(i-1) + three * u(i  ))**2
      
      alpha0 = d0 / (epsilon + beta0)**2
      alpha1 = d1 / (epsilon + beta1)**2
      alpha2 = d2 / (epsilon + beta2)**2
      
      alphasum = alpha0 + alpha1 + alpha2
      
      w0 = alpha0 / alphasum
      w1 = alpha1 / alphasum
      w2 = alpha2 / alphasum
      
      u0plushalf = c00 * u(i  ) + c01 * u(i+1) + c02 * u(i+2)
      u1plushalf = c10 * u(i-1) + c11 * u(i  ) + c12 * u(i+1)
      u2plushalf = c20 * u(i-2) + c21 * u(i-1) + c22 * u(i  )
      
      uminus(i) = w0 * u0plushalf + &
           w1 * u1plushalf + &
           w2 * u2plushalf
      
    end do
    
  end subroutine WENO3_Left
  
  subroutine WENO3_Right(np, ngz, u, uplus)
    
    use real_type_mod
    
    implicit none
    
    integer :: i, np, ngz
    real(kind=wp), dimension(1 - ngz : np + ngz) :: u, uplus
    real(kind=wp) :: one, two, three, four, five, six, seven, &
         ten, eleven, twelve, thirteen, &
         ThirteenByTwelve, Quarter
    parameter (one = 1)
    parameter (two = 2)
    parameter (three = 3)
    parameter (four = 4)
    parameter (five = 5)
    parameter (six = 6)
    parameter (seven = 7)
    parameter (ten = 10)
    parameter (eleven = 11)
    parameter (twelve = 12)
    parameter (thirteen = 13)
    parameter (ThirteenByTwelve = thirteen / twelve)
    parameter (Quarter = one / four)
    real(kind=wp) :: dtilde0, dtilde1, dtilde2
    parameter (dtilde0 = one / ten)
    parameter (dtilde1 = six / ten)
    parameter (dtilde2 = three / ten)
    real(kind=wp) :: ctilde00,ctilde01,ctilde02,ctilde10,ctilde11,&
         ctilde12,ctilde20,ctilde21,ctilde22
    parameter (ctilde00 = eleven / six)
    parameter (ctilde01 = -seven / six)
    parameter (ctilde02 = two / six)
    parameter (ctilde10 = two / six)
    parameter (ctilde11 = five / six)
    parameter (ctilde12 = -one / six)
    parameter (ctilde20 = -one / six)
    parameter (ctilde21 = five / six)
    parameter (ctilde22 = two / six)
    
    real(kind=wp) :: betatilde0, betatilde1, betatilde2
    real(kind=wp) :: epsilon
    real(kind=wp) :: alphatilde0, alphatilde1, alphatilde2, alphatildesum
    real(kind=wp) :: wtilde0, wtilde1, wtilde2
    real(kind=wp) :: u0minushalf, u1minushalf, u2minushalf
    
    epsilon = 1.d-6
    
    do i = 0, np+1
      
      betatilde0 = ThirteenByTwelve * (u(i  ) - two * u(i+1) + u(i+2))**2 + &
           Quarter * (three * u(i  ) - four * u(i+1) +         u(i+2))**2
      betatilde1 = ThirteenByTwelve * (u(i-1) - two * u(i  ) + u(i+1))**2 + &
           Quarter * (        u(i-1)                 -         u(i+1))**2
      betatilde2 = ThirteenByTwelve * (u(i-2) - two * u(i-1) + u(i  ))**2 + &
           Quarter * (        u(i-2) - four * u(i-1) + three * u(i  ))**2
      
      alphatilde0 = dtilde0 / (epsilon + betatilde0)**2
      alphatilde1 = dtilde1 / (epsilon + betatilde1)**2
      alphatilde2 = dtilde2 / (epsilon + betatilde2)**2
      
      alphatildesum = alphatilde0 + alphatilde1 + alphatilde2
      
      wtilde0 = alphatilde0 / alphatildesum
      wtilde1 = alphatilde1 / alphatildesum
      wtilde2 = alphatilde2 / alphatildesum
      
      u0minushalf = ctilde00 * u(i  ) + ctilde01 * u(i+1) + ctilde02 * u(i+2)
      u1minushalf = ctilde10 * u(i-1) + ctilde11 * u(i  ) + ctilde12 * u(i+1)
      u2minushalf = ctilde20 * u(i-2) + ctilde21 * u(i-1) + ctilde22 * u(i  )
      
      uplus(i) = wtilde0 * u0minushalf + &
           wtilde1 * u1minushalf + &
           wtilde2 * u2minushalf
      
    end do
    
  end subroutine WENO3_Right
  
end module reconstruct_mod
