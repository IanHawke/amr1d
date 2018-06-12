!!$The computation of the flux; effectively the Riemann solver

subroutine flux_formula(nv, np, ngz, ul, ur, flux)

  use real_type_mod
  use system_parameters_mod
  use method_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: ul, ur
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(out) :: flux

!!$  The initial states for the Riemann problem; avoids confusion

  real(kind=wp), dimension(nv, 1-ngz:np+ngz) :: ql, qr

!!$  The left and right fluxes

  real(kind=wp), dimension(nv, 1-ngz:np+ngz) :: fl, fr

  integer :: i, k

!!$  Temporaries for the characteristic speed calculation

  real(kind=wp), dimension(1-ngz:np+ngz) :: char_l, char_r, char_max

!!$  Godunov Riemann problem solution

  flux = 0.0_wp

!!$  Advection equation

  if (advection) then

    if (a > 0.0_wp) then
      
      flux(:, 1:np+1) = a * ur(:, 0:np)
      
    else
      
      flux(:, 1:np+1) = a * ul(:, 1:np+1)
      
    end if

!!$  Burgers equation

  else if (burgers) then

    ql(:, 1:np+1) = ur(:, 0:np)
    qr(:, 1:np+1) = ul(:, 1:np+1)

    do i = 1, np + 1
      
      do k = 1, nv

        if (ql(k, i) < qr(k, i)) then
        
          if (ql(k, i) > 0.0_wp) then
            
            flux(k, i) = 0.5_wp * ql(k, i)**2
            
          else if (ql(k, i) * qr(k, i) < 0.0_wp) then
            
            flux(k, i) = 0.0_wp
            
          else
            
            flux(k, i) = 0.5_wp * qr(k, i)**2
            
          end if
          
        else
          
          if (ql(k, i) + qr(k, i) > 0.0_wp) then
            
            flux(k, i) = 0.5_wp * ql(k, i)**2
            
          else
            
            flux(k, i) = 0.5_wp * qr(k, i)**2
            
          end if
          
        end if
      
      end do
      
    end do

  else if (scalarwave) then

    ql(:, 1:np+1) = ur(:, 0:np)
    qr(:, 1:np+1) = ul(:, 1:np+1)

    call pointwise_flux(nv, np, ngz, ql, fl)
    call pointwise_flux(nv, np, ngz, qr, fr)

!!$    Use an HLLE approximate Riemann Solver

    flux = 0.0_wp

    char_max = 1.0_wp

    call hlle(nv, nv, np, ngz, &
           ql, qr, &
           fl, fr, &
           char_max, &
           flux)
      
  else if (euler) then

    ql(:, 1:np+1) = ur(:, 0:np)
    qr(:, 1:np+1) = ul(:, 1:np+1)

    call pointwise_flux(nv, np, ngz, ql, fl)
    call pointwise_flux(nv, np, ngz, qr, fr)

!!$    Use an HLLE approximate Riemann Solver

    flux = 0.0_wp

!!$    Finding a robust maximum characteristic speed is not easy
!!$    in extreme cases. This should work for TVD at least...
!!$    Essentially what is required is the maximum over the
!!$    _relevant_ range of the values. This may stretch over
!!$    at least two cells

    char_max = 1.0_wp

    call max_char(nv, np, ngz, ul, char_l)
    call max_char(nv, np, ngz, ur, char_r)

    do i = 1, np + 1
      char_max(i) = max(maxval(char_l(i:i+1)), &
                        maxval(char_r(i-1:i)))
    end do

    if (rs_HLLE) then

!!$      Only apply the flux formula to the conserved variables

      call hlle(nv, 3, np, ngz, &
           ql, qr, &
           fl, fr, &
           char_max, &
           flux)

!!$    HLLC; here it requires the primitive variables as well.

    else if (rs_HLLC) then

      call hllc(nv, 6, np, ngz, &
           ql, qr, &
           fl, fr, &
           char_max, &
           flux)

    else

      write(*,*) 'Do not understand which Riemann solver you want'
      STOP

    end if
    
!!$    write(*,*) 'flux', flux(3, :)
!!$    STOP
  end if
  
end subroutine flux_formula

!!$The pointwise flux, useful in various places

subroutine pointwise_flux(nv, np, ngz, u, flux)

  use real_type_mod
  use system_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: u
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(out) :: flux

  flux = 0.0_wp

  if (advection) then

    flux = a * u

  else if (burgers) then
    
    flux = 0.5_wp * u**2

  else if (scalarwave) then

    flux(1, :) = 0.0_wp
    flux(2, :) = - u(3, :)
    flux(3, :) = - u(2, :)

  else if (euler) then

    flux(1, :) = u(2, :)
    flux(2, :) = u(2, :) * u(5, :) + u(6, :)
    flux(3, :) = (u(3, :) + u(6, :)) * u(5, :)

  end if
  
end subroutine pointwise_flux

!!$HLLE flux formula

subroutine hlle(nv, nc, np, ngz, ql, qr, fl, fr, char_max, flux)

  use real_type_mod
  implicit none

  integer, intent(in) :: nv, nc, np, ngz
  real(kind=wp), intent(in), dimension(nv, 1-ngz:np+ngz) :: ql, qr, fl, fr
  real(kind=wp), intent(in), dimension(1-ngz:np+ngz) :: char_max
  real(kind=wp), intent(out), dimension(nv, 1-ngz:np+ngz) :: flux

  integer :: k

  do k = 1, nc

    flux(k, 1:np+1) =       0.5_wp * (fl(k, 1:np+1) + fr(k, 1:np+1) + &
                  char_max(1:np+1) * (ql(k, 1:np+1) - qr(k, 1:np+1)))

  end do
  
end subroutine hlle

!!$HLLC flux formula - Euler equations only

subroutine hllc(nv, nc, np, ngz, ql, qr, fl, fr, char_max, flux)

  use real_type_mod
  implicit none

  integer, intent(in) :: nv, nc, np, ngz
  real(kind=wp), intent(in), dimension(nv, 1-ngz:np+ngz) :: ql, qr, fl, fr
  real(kind=wp), intent(in), dimension(1-ngz:np+ngz) :: char_max
  real(kind=wp), intent(out), dimension(nv, 1-ngz:np+ngz) :: flux

  real(kind=wp) :: star_char
  real(kind=wp), dimension(nv) :: qstarl, qstarr

  integer :: k, i

  do i = 1, np + 1

!!$      See eq 10.58 of Toro

    star_char = ( qr(6, i) - ql(6, i) + &
         ql(4, i) * ql(5, i) * (-char_max(i) - ql(5, i)) - &
         qr(4, i) * qr(5, i) * ( char_max(i) - qr(5, i)) ) / &
         (ql(4, i) * (-char_max(i) - ql(5, i)) - &
          qr(4, i) * ( char_max(i) - qr(5, i)))

    qstarl = 0.0_wp
    qstarr = 0.0_wp

    qstarl(1) = 1.0_wp
    qstarl(2) = star_char
    qstarl(3) = ql(3,i) / ql(1,i) + (star_char - ql(5,i)) * &
         (star_char + ql(6,i) / ql(1,i) / (-char_max(i) - ql(5,i)))
    qstarl = qstarl * ql(1,i) * (-char_max(i) - ql(5,i)) / (-char_max(i) - star_char)

    qstarr(1) = 1.0_wp
    qstarr(2) = star_char
    qstarr(3) = qr(3,i) / qr(1,i) + (star_char - qr(5,i)) * &
         (star_char + qr(6,i) / qr(1,i) / ( char_max(i) - qr(5,i)))
    qstarr = qstarr * qr(1,i) * ( char_max(i) - qr(5,i)) / ( char_max(i) - star_char)

    do k = 1, nc

      if (star_char > 0.0_wp) then

        flux(k, i) = fl(k, i) - char_max(i) * (qstarl(k) - ql(k, i))

      else

        flux(k, i) = fr(k, i) + char_max(i) * (qstarr(k) - qr(k, i))

      end if
      
    end do
    
  end do
  
end subroutine hllc

!!$Find maximum characteristic speed

subroutine max_char(nv, np, ngz, u, char_max)

  use real_type_mod
  use system_parameters_mod
  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), intent(in), dimension(nv, 1-ngz:np+ngz) :: u
  real(kind=wp), intent(out), dimension(1-ngz:np+ngz) :: char_max

  real(kind=wp), dimension(1-ngz:np+ngz) :: cs

  char_max = 1.0_wp

  if (advection) then

    char_max = a

  else if (burgers) then

    char_max = maxval(abs(u))

  else if (scalarwave) then

!!$    Fall through as char_max = 1

  else if (euler) then

    cs(1:np+1) = sqrt(euler_gamma * u(6, 1:np+1) / u(4, 1:np+1))

    char_max(1:np+1) = abs(u(5, 1:np+1)) + cs(1:np+1)

  end if

end subroutine max_char


subroutine maxval_char(nv, np, ngz, u, char_max)

  use real_type_mod
  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), intent(in), dimension(nv, 1-ngz:np+ngz) :: u
  real(kind=wp), intent(out) :: char_max

  real(kind=wp), dimension(1-ngz:np+ngz) :: char

  call max_char(nv, np, ngz, u, char)

  char_max = maxval(char)

end subroutine maxval_char

subroutine eigenvalues(nv, np, ngz, u, evals)

  use real_type_mod
  use system_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: u
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(out) :: evals

  real(kind=wp), dimension(:), allocatable :: cs

  integer :: ierr

  evals = 1.0_wp

  if (advection) then

    evals = a
    
  else if (burgers) then

    evals = u

  else if (scalarwave) then

    evals(1, :) =  0.0_wp
    evals(2, :) = -1.0_wp
    evals(3, :) =  1.0_wp

  else if (euler) then

    call conservative2primitive(nv, np, ngz, u)

    allocate(cs   (    1-ngz:np+ngz), &
             STAT = ierr)
    if (ierr .ne. 0) then
      write(*,*) 'Error allocating errors in eigenvectors!'
      STOP
    end if
        
    evals(2, :) = u(5, :)
    cs(:) = sqrt(euler_gamma * u(6, :) / u(4, :))
    evals(1, :) = evals(2, :) - cs(:)
    evals(3, :) = evals(2, :) + cs(:)

    deallocate(cs)


  else

    write(*,*) 'What system are you evolving? eigenvalues is confused'
    STOP

  end if
  
end subroutine eigenvalues

subroutine eigenvectors(nv, np, ngz, u, leftev, rightev)

  use real_type_mod
  use system_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: u
  real(kind=wp), dimension(nv, nv, 1-ngz:np+ngz), intent(out) :: leftev, rightev
  
  real(kind=wp), dimension(:, :), allocatable :: evals
  real(kind=wp), dimension(:), allocatable :: cs, h, b1, b2
  real(kind=wp) :: gm1
  
  integer :: ierr, i, var
  
  leftev  = 0.0_wp
  rightev = 0.0_wp
  
  if (euler) then
  
    gm1 = euler_gamma - 1.0_wp
  
    allocate(evals(nv, 1-ngz:np+ngz), &
             cs   (    1-ngz:np+ngz), &
             h    (    1-ngz:np+ngz), &
             b1   (    1-ngz:np+ngz), &
             b2   (    1-ngz:np+ngz), &
             STAT = ierr)
    if (ierr .ne. 0) then
      write(*,*) 'Error allocating errors in eigenvectors!'
      STOP
    end if
    
    evals(2, :) = u(5, :)
    cs(:) = sqrt(euler_gamma * u(6, :) / u(4, :))
    evals(1, :) = evals(2, :) - cs(:)
    evals(3, :) = evals(2, :) + cs(:)
    
    h(:) = cs(:)**2 / gm1 + 0.5_wp * u(5, :)**2
    b1(:) = gm1 / cs(:)**2
    b2(:) = 0.5_wp * b1 * u(5, :)**2
    
    rightev(1, 1, :) = 1.0_wp
    rightev(1, 2, :) = evals(1, :)
    rightev(1, 3, :) = h(:) - cs(:) * u(5, :)
    rightev(2, 1, :) = 1.0_wp
    rightev(2, 2, :) = u(5, :)
    rightev(2, 3, :) = 0.5_wp * u(5, :) * u(5, :)
    rightev(3, 1, :) = 1.0_wp
    rightev(3, 2, :) = evals(3, :)
    rightev(3, 3, :) = h(:) + cs(:) * u(5, :)
    
    leftev(1, 1, :) = 0.5_wp * (b2(:) + u(5, :) / cs(:))
    leftev(1, 2, :) = -0.5_wp * (b1(:) * u(5, :) + 1.0_wp / cs(:))
    leftev(1, 3, :) = 0.5_wp * b1(:)
    leftev(2, 1, :) = 1.0_wp - b2(:)
    leftev(2, 2, :) = b1(:) * u(5, :)
    leftev(2, 3, :) = -b1(:)
    leftev(3, 1, :) = 0.5_wp * (b2(:) - u(5, :) / cs(:))     
    leftev(3, 2, :) = 0.5_wp * (-b1(:) * u(5, :) + 1.0_wp / cs(:))
    leftev(3, 3, :) = 0.5_wp * b1(:)

    deallocate(evals, &
               cs, &
               h , &
               b1, &
               b2)
    
  else

!!$    At present only the Euler system has non-trivial values

    do var = 1, nv
      leftev (var, var, :) = 1.0_wp
      rightev(var, var, :) = 1.0_wp
    end do
    
  end if

!!$  Leave in debugging code to check that the eigenvectors are
!!$  orthonormal

!!$    do i = 1, np
!!$      write(*,*) i, sum(leftev(1, :, i)*rightev(1, :, i)), &
!!$                    sum(leftev(2, :, i)*rightev(1, :, i)), &
!!$                    sum(leftev(3, :, i)*rightev(1, :, i)), &
!!$                    sum(leftev(1, :, i)*rightev(2, :, i)), &
!!$                    sum(leftev(2, :, i)*rightev(2, :, i)), &
!!$                    sum(leftev(3, :, i)*rightev(2, :, i)), &
!!$                    sum(leftev(1, :, i)*rightev(3, :, i)), &
!!$                    sum(leftev(2, :, i)*rightev(3, :, i)), &
!!$                    sum(leftev(3, :, i)*rightev(3, :, i))
!!$    end do
!!$    STOP
  
end subroutine eigenvectors
