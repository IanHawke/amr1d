module sources_mod

  use real_type_mod
  use system_parameters_mod

  contains

!!$Source terms

    subroutine sources(ngz, u, source)

      implicit none

      integer, intent(in) :: ngz
      real(kind=wp), dimension(:, 1-ngz:), intent(in) :: u
      real(kind=wp), dimension(:, 1-ngz:), intent(out) :: source

      source = 0.0_wp

      if (scalarwave) then
        
        source(1, :) = u(2, :)

      end if

    end subroutine sources

end module
