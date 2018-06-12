!!$Output routine

subroutine output(name, nv, np, ngz, x, u, u_exact, error)

  use real_type_mod
  use grid_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(1-ngz:np+ngz), intent(in) :: &
       x
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: &
       u, u_exact, error

  character(len=*), intent(in) :: name
  character(len=200) :: filestring, filestring_number, filestring_number_trim
  integer :: filestring_len, i, k

  write(filestring_number, *) np
  filestring_len = int(log10(dble(np)))+1
  if (filestring_len > 9) then
    write(*,*) 'Too many points for file output!'
    stop
  end if
  filestring_number_trim = trim( &
       filestring_number(len_trim(filestring_number) + 1 - filestring_len:&
       len_trim(filestring_number)))
  filestring = trim(name//filestring_number_trim(1:filestring_len)//'.dat')

  open(10, file=filestring, status='replace')

  if (output_exact_and_errors) then

    do i = 1, np
      write(10, '(es21.12E3,a1)',advance='NO') x(i), ' '
      do k = 1, nv
        write(10, '(3es21.12E3,a1)',advance='NO') u(k,i), &
             u_exact(k,i), error(k,i)
      end do
      write(10, *)
    end do

  else

    do i = 1, np
      write(10, '(es21.12E3,a1)',advance='NO') x(i), ' '
      do k = 1, nv
        write(10, '(es21.12E3,a1)',advance='NO') u(k,i), ' '
      end do
      write(10, *)
    end do

  end if

  close(10)
  
end subroutine output

subroutine output_errors(name, nv, np, ngz, n, dx, error)

  use real_type_mod
  use grid_parameters_mod

  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), intent(in) :: dx
  real(kind=wp), dimension(nv, 1-ngz:np+ngz), intent(in) :: error

  character(len=*), intent(in) :: name
  character(len=200) :: filestring
  integer :: n

  filestring = trim(name//'.errors')

  if (n == 1) then
    open(10, file=filestring, status='replace')
    write(10, '(a23)') '#  dx  L_1  L_2  L_inf'
  else
    open(10, file=filestring, status='old', access='append')
  end if
  
  write(10, '(4es21.12E3)') dx, &
       dx * sum(error), &
       sqrt( dx * sum(error**2) ), &
       maxval(error)

  close(10)

end subroutine output_errors
