!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Adjoint of poly1d_vert_w3_reconstruction_kernel. Added to remove necessity
!!        of global variable "stencil" (which exists in tangent linear kernel).
!!        Fixed by PSyclone issue #1312.
module invoke_atl_poly1d_vert_w3_recon_mod

  use field_mod,                    only : field_type, field_proxy_type
  use r_tran_field_mod,             only : r_tran_field_type, r_tran_field_proxy_type
  use mesh_mod,                     only : mesh_type
  use constants_mod,                only : r_def, i_def, r_double, r_solver, r_tran, l_def
  use r_solver_field_mod,           only : r_solver_field_type, &
                                           r_solver_field_proxy_type

  implicit none
  public :: invoke_atl_poly1d_vert_w3_recon_kernel_type

contains


    SUBROUTINE invoke_atl_poly1d_vert_w3_recon_kernel_type(field_new, field_old, ls_field_old, vert_flux_coeffs, ndata_v, &
                                                           vertical_order, logspace, stencil)
      USE atl_poly1d_vert_w3_reconstruction_kernel_mod, ONLY: atl_poly1d_vert_w3_reconstruction_code
      USE mesh_mod, ONLY: mesh_type
      INTEGER(KIND=i_def), intent(in) :: ndata_v, vertical_order
      LOGICAL(KIND=l_def), intent(in) :: logspace
      TYPE(field_type), intent(in) :: field_new, field_old, ls_field_old, vert_flux_coeffs
      integer(kind=i_def), intent(inout), allocatable, dimension(:,:,:) :: stencil
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers_field_new
      REAL(KIND=r_def), pointer, dimension(:) :: vert_flux_coeffs_data
      REAL(KIND=r_def), pointer, dimension(:) :: ls_field_old_data
      REAL(KIND=r_def), pointer, dimension(:) :: field_old_data
      REAL(KIND=r_def), pointer, dimension(:) :: field_new_data
      TYPE(field_proxy_type) field_new_proxy, field_old_proxy, ls_field_old_proxy, vert_flux_coeffs_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_field_new(:,:), map_adspc2_vert_flux_coeffs(:,:), &
&map_w3(:,:)
      INTEGER(KIND=i_def) ndf_adspc1_field_new, undf_adspc1_field_new, ndf_w3, undf_w3, ndf_adspc2_vert_flux_coeffs, &
&undf_adspc2_vert_flux_coeffs
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh
      !
      ! Initialise field and/or operator proxies
      !
      field_new_proxy = field_new%get_proxy()
      field_new_data => field_new_proxy%data
      field_old_proxy = field_old%get_proxy()
      field_old_data => field_old_proxy%data
      ls_field_old_proxy = ls_field_old%get_proxy()
      ls_field_old_data => ls_field_old_proxy%data
      vert_flux_coeffs_proxy = vert_flux_coeffs%get_proxy()
      vert_flux_coeffs_data => vert_flux_coeffs_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers_field_new = field_new_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => field_new_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_field_new => field_new_proxy%vspace%get_whole_dofmap()
      map_w3 => field_old_proxy%vspace%get_whole_dofmap()
      map_adspc2_vert_flux_coeffs => vert_flux_coeffs_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_field_new
      !
      ndf_adspc1_field_new = field_new_proxy%vspace%get_ndf()
      undf_adspc1_field_new = field_new_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for w3
      !
      ndf_w3 = field_old_proxy%vspace%get_ndf()
      undf_w3 = field_old_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_vert_flux_coeffs
      !
      ndf_adspc2_vert_flux_coeffs = vert_flux_coeffs_proxy%vspace%get_ndf()
      undf_adspc2_vert_flux_coeffs = vert_flux_coeffs_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(cell)
      !$omp do schedule(static)
      DO cell = loop0_start, loop0_stop, 1
        CALL atl_poly1d_vert_w3_reconstruction_code(nlayers_field_new, field_new_data, field_old_data, ls_field_old_data, &
&vert_flux_coeffs_data, ndata_v, vertical_order, logspace, ndf_adspc1_field_new, undf_adspc1_field_new, &
&map_adspc1_field_new(:,cell), ndf_w3, undf_w3, map_w3(:,cell), ndf_adspc2_vert_flux_coeffs, undf_adspc2_vert_flux_coeffs, &
&map_adspc2_vert_flux_coeffs(:,cell), stencil)
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL field_new_proxy%set_dirty()
      CALL field_old_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_atl_poly1d_vert_w3_recon_kernel_type

end module invoke_atl_poly1d_vert_w3_recon_mod
