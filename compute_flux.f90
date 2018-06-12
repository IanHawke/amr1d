!!$Compute the flux

subroutine intercell_flux(nv, np, ngz, u, grid_data)

  use real_type_mod
  use data_mod
  use method_parameters_mod
  use reconstruct_mod
  use flux_splitting_mod
  
  implicit none

  integer, intent(in) :: nv, np, ngz
  real(kind=wp), dimension(nv, 1 - ngz:np + ngz), intent(in) :: u
  type(data) :: grid_data

  real(kind=wp), dimension(:, :), allocatable :: tmp_flux
  integer :: ierr

  if (Finite_Volume) then

!!$  Reconstruct to cell boundaries

    call reconstruct(nv, np, ngz, u, grid_data)

!!$  Solve the Riemann problem

    call flux_formula(nv, np, ngz, &
         grid_data%ul, grid_data%ur, grid_data%flux)

  else if (Finite_Difference) then

    call compute_flux_weno_fluxsplit(nv, np, ngz, u, grid_data)

!!$    Convert to characteristic variables

!!$    We actually have to do this on an interface by interface
!!$    basis, so the following is all wrong...

!!$    allocate(tmp_flux(nv, 1-ngz:np+ngz), &
!!$             STAT=ierr)
!!$    if (ierr .ne. 0) then
!!$      write(*,*) 'Failed to allocate array in intercell_flux!'
!!$      STOP
!!$    end if
!!$    
!!$    call pointwise_flux(nv, np, ngz, u, tmp_flux)
!!$
!!$    call convert_to_characteristics(nv, np, ngz, &
!!$         u, tmp_flux, grid_data)
!!$
!!$    deallocate(tmp_flux)

!!$    Do the flux splitting

!!$    call set_flux_splitting(nv, np, ngz, &
!!$         u, grid_data)

!!$    Reconstruct the split fluxes

!!$    call reconstruct_split_flux(nv, np, ngz, grid_data)

!!$    Convert back to physical space

!!$    call convert_from_characteristics(nv, np, ngz, &
!!$         u, grid_data)

!!$    Recombine to get the final result

!!$    call recombine_flux_splitting(nv, np, ngz, &
!!$         grid_data)

  else

    write(*,*) "Finite volume or difference?"
    write(*,*) "intercell_flux is confused."
    STOP
  end if
  
end subroutine intercell_flux
