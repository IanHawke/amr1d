!!$Compute a flux splitting and later recombine it

module flux_splitting_mod

  use real_type_mod
  use data_mod
  use method_parameters_mod

contains

  subroutine set_flux_splitting(nv, np, ngz, u, f, f_plus, f_minus)

    implicit none

    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, -2:3, 1 - ngz:np + ngz), intent(in) :: u, f
    real(kind=wp), dimension(nv, -2:3, 1 - ngz:np + ngz), intent(out) :: f_plus, &
                                                                         f_minus

    real(kind=wp), dimension(:, :), allocatable :: evals
    integer :: ierr, i, var, imin, imax

    real(kind=wp) :: alpha
    
!!$    call maxval_char(nv, np, ngz, u, alpha)
!!$
!!$    f_minus = 0.5_wp * (f - alpha * u)
!!$    f_plus  = 0.5_wp * (f + alpha * u)

    allocate(evals(nv, 1 - ngz:np + ngz), &
             STAT = ierr)
    if (ierr .ne. 0) then
      write(*,*) 'Failed to allocate memory in set_flux_splitting'
      STOP
    end if
    
    call eigenvalues(nv, np, ngz, u, evals)

    do i = 1 - ngz, np + ngz
      do var = 1, nv
        imin = max(i - 3, 1 - ngz)
        imax = min(i + 3, np + ngz)

!!$        We could take the maximum over just this characteristic field.
!!$        However, that seems to lead to overshoots near tails of rarefactions.
!!$        There may be an easy fix for that, but it is easiest to take the 
!!$        maximum over all characteristic fields.

        alpha = maxval(abs(evals(:, imin:imax)))

        f_minus(var, :, i) = 0.5_wp * (f(var, :, i) - &
                               alpha * u(var, :, i))
        f_plus (var, :, i) = 0.5_wp * (f(var, :, i) + &
                               alpha * u(var, :, i))
      end do
    end do
    
    deallocate(evals)

  end subroutine set_flux_splitting
  
  subroutine recombine_flux_splitting(nv, np, ngz, &
                                      f_plus_right, f_minus_left, f)
  
    implicit none

    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(in) :: f_plus_right, &
                                                                  f_minus_left
    real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(out) :: f

    f(:, 1:np + 1) = &
         f_minus_left(:, 1:np + 1) + &
         f_plus_right(:, 0:np    )

  end subroutine recombine_flux_splitting

  subroutine set_eigenvectors(nv, np, ngz, u, grid_data)

    implicit none

    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(in) :: u
    type(data) :: grid_data
    
    real(kind=wp), dimension(:, :), allocatable :: uint
    integer :: ierr, var, i, k, l, lmin, lmax

    allocate(uint   (nv, 1-ngz:np+ngz), &
             STAT = ierr)
    if (ierr .ne. 0) then
      write(*,*) 'Allocation problems in convert_to_characteristics'
      STOP
    end if

!!$    Do not set to zero or there will be zero density problems...

    uint = 1.0_wp

    do i = 1 - ngz, np + ngz - 1
!!$      Simple arithmetic average for the interface state
      uint(:, i) = 0.5_wp * (u(:, i) + u(:, i+1))
    end do

    call conservative2primitive(nv, np, ngz, uint)
    
    call eigenvectors(nv, np, ngz, &
         uint, grid_data%leftev, grid_data%rightev)

    deallocate(uint)

  end subroutine set_eigenvectors
  
  subroutine convert_to_characteristics(nv, np, ngz, u, f, leftev, u_char, f_char)

    implicit none

    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(in) :: u, f
    real(kind=wp), dimension(nv, nv, 1 - ngz:np + ngz), intent(in) :: leftev
    real(kind=wp), dimension(nv, -2:3, 1 - ngz:np + ngz), intent(out) :: u_char, f_char

    integer :: var, i, k, l, lmin, lmax

    u_char = 0.0_wp
    f_char = 0.0_wp

    do i = 1 - ngz, np + ngz
      lmin = max(1  - ngz, i - 2) - i
      lmax = min(np + ngz, i + 3) - i
      do l = lmin, lmax
        do var = 1, nv
          do k = 1, nv
            u_char(var, l, i) = u_char(var, l, i) + &
                 leftev(var, k, i) * u(k, i + l)
            f_char(var, l, i) = f_char(var, l, i) + &
                 leftev(var, k, i) * f(k, i + l)
          end do
        end do
      end do
      do var = 1, nv
        u_char(var, -2:lmin - 1, i) = u_char(var, lmin, i)
        f_char(var, -2:lmin - 1, i) = f_char(var, lmin, i)
        u_char(var, lmax + 1:3 , i) = u_char(var, lmax, i)
        f_char(var, lmax + 1:3 , i) = f_char(var, lmax, i)
      end do
    end do

  end subroutine convert_to_characteristics
  
  subroutine convert_from_characteristics(nv, np, ngz, &
                                          f_char_plus_right, f_char_minus_left, &
                                          rightev, &
                                          f_plus_right, f_minus_left)

    implicit none

    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(in) :: f_char_plus_right, &
                                                                  f_char_minus_left
    real(kind=wp), dimension(nv, nv, 1 - ngz:np + ngz), intent(in) :: rightev
    real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(out) :: f_plus_right, &
                                                                   f_minus_left

    integer :: var, i, k

    do i = 1 - ngz, np + ngz - 1
      do var = 1, nv
        f_minus_left(var, i + 1) = 0.0_wp
        f_plus_right(var, i    ) = 0.0_wp
        do k = 1, nv
          f_minus_left(var, i + 1) = f_minus_left(var, i + 1) + &
               rightev(k, var, i) * f_char_minus_left(k, i + 1)
          f_plus_right(var, i    ) = f_plus_right(var, i    ) + &
               rightev(k, var, i) * f_char_plus_right(k, i    )
        end do
      end do
    end do

  end subroutine convert_from_characteristics

  subroutine compute_flux_weno_fluxsplit(nv, np, ngz, u, grid_data)
  
    implicit none

    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(in) :: u
    type(data) :: grid_data

!!$    Firstly we compute the flux at the grid points.
!!$    Use the flux array as a temporary.

    call pointwise_flux(nv, np, ngz, u, grid_data%flux)

!!$    Set up the eigenvectors so that characteristic splitting can
!!$    be performed

    call set_eigenvectors(nv, np, ngz, u, grid_data)

!!$    Convert to characteristic variables

    call convert_to_characteristics(nv, np, ngz, &
         u, grid_data%flux, &
         grid_data%leftev, &
         grid_data%u_char, grid_data%f_char)

!!$    Split the flux using the LF flux splitting

    call set_flux_splitting(nv, np, ngz, &
         grid_data%u_char, grid_data%f_char, &
         grid_data%f_char_plus, grid_data%f_char_minus)

!!$    Use the WENO algorithm to reconstruct at each interface

    call weno_fs_reconstruct(nv, np, ngz, &
         grid_data%f_char_plus, grid_data%f_char_minus, &
         grid_data%f_char_plus_right, grid_data%f_char_minus_left)

!!$    Convert back to physical space

    call convert_from_characteristics(nv, np, ngz, &
         grid_data%f_char_plus_right, grid_data%f_char_minus_left, &
         grid_data%rightev, &
         grid_data%f_plus_right, grid_data%f_minus_left)

!!$    Recombine the split fluxes to get the final answer

    call recombine_flux_splitting(nv, np, ngz, &
         grid_data%f_plus_right, grid_data%f_minus_left, &
         grid_data%flux)

  end subroutine compute_flux_weno_fluxsplit

  subroutine weno_fs_reconstruct(nv, np, ngz, &
       f_plus, f_minus, &
       f_plus_right, f_minus_left)

    implicit none

    integer, intent(in) :: nv, np, ngz
    real(kind=wp), dimension(nv, -2:3, 1 - ngz:np + ngz), intent(in) :: f_plus, f_minus
    real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(out) :: f_plus_right, &
                                                                   f_minus_left

    integer :: var, i

    do i = 0, np
      do var = 1, nv
        call WENO5_Left_point (f_plus (var, -2:2, i), f_plus_right(var, i    ))
        call WENO5_Right_point(f_minus(var, -1:3, i), f_minus_left(var, i + 1))
      end do
    end do
    
  end subroutine weno_fs_reconstruct
  
  subroutine WENO5_Left_point(u, uminus)
    
    use real_type_mod
    
    implicit none
    
    real(kind=wp), dimension(5), intent(in) :: u
    real(kind=wp), intent(out) :: uminus
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
    
    beta0 = ThirteenByTwelve * (u(3  ) - two * u(4) + u(5))**2 + &
         Quarter * (three * u(3) - four * u(4) +         u(5))**2
    beta1 = ThirteenByTwelve * (u(2) - two * u(3) + u(4))**2 + &
         Quarter * (        u(2)                 -       u(4))**2
    beta2 = ThirteenByTwelve * (u(1) - two * u(2) + u(3))**2 + &
         Quarter * (        u(1) - four * u(2) + three * u(3))**2
      
    alpha0 = d0 / (epsilon + beta0)**2
    alpha1 = d1 / (epsilon + beta1)**2
    alpha2 = d2 / (epsilon + beta2)**2
      
    alphasum = alpha0 + alpha1 + alpha2
      
    w0 = alpha0 / alphasum
    w1 = alpha1 / alphasum
    w2 = alpha2 / alphasum
      
    u0plushalf = c00 * u(3) + c01 * u(4) + c02 * u(5)
    u1plushalf = c10 * u(2) + c11 * u(3) + c12 * u(4)
    u2plushalf = c20 * u(1) + c21 * u(2) + c22 * u(3)
      
    uminus = w0 * u0plushalf + &
         w1 * u1plushalf + &
         w2 * u2plushalf
          
  end subroutine WENO5_Left_point
  
  subroutine WENO5_Right_point(u, uplus)
    
    use real_type_mod
    
    implicit none
    
    real(kind=wp), dimension(5), intent(in) :: u
    real(kind=wp), intent(out) :: uplus
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
    
    betatilde0 = ThirteenByTwelve * (u(3) - two * u(4) + u(5))**2 + &
         Quarter * (three * u(3) - four * u(4) +         u(5))**2
    betatilde1 = ThirteenByTwelve * (u(2) - two * u(3) + u(4))**2 + &
         Quarter * (        u(2)                 -         u(4))**2
    betatilde2 = ThirteenByTwelve * (u(1) - two * u(2) + u(3))**2 + &
         Quarter * (        u(1) - four * u(2) + three * u(3))**2
    
    alphatilde0 = dtilde0 / (epsilon + betatilde0)**2
    alphatilde1 = dtilde1 / (epsilon + betatilde1)**2
    alphatilde2 = dtilde2 / (epsilon + betatilde2)**2
    
    alphatildesum = alphatilde0 + alphatilde1 + alphatilde2
    
    wtilde0 = alphatilde0 / alphatildesum
    wtilde1 = alphatilde1 / alphatildesum
    wtilde2 = alphatilde2 / alphatildesum
    
    u0minushalf = ctilde00 * u(3) + ctilde01 * u(4) + ctilde02 * u(5)
    u1minushalf = ctilde10 * u(2) + ctilde11 * u(3) + ctilde12 * u(4)
    u2minushalf = ctilde20 * u(1) + ctilde21 * u(2) + ctilde22 * u(3)
    
    uplus = wtilde0 * u0minushalf + &
         wtilde1 * u1minushalf + &
         wtilde2 * u2minushalf
    
  end subroutine WENO5_Right_point
  
end module flux_splitting_mod
