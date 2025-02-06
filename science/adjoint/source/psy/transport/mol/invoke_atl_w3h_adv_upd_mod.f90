!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Regular invocation of the line-by-line kernel has incorrect
!!        halo depths, so need PSyKAl lite solution. Addressed by
!!        PSyclone issue #2886.
MODULE invoke_atl_w3h_adv_upd_mod

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: invoke_atl_w3h_adv_upd_kernel_type

CONTAINS

    SUBROUTINE invoke_atl_w3h_adv_upd_kernel_type( advective_increment, tracer, tracer_stencil_extent, &
                                                   ls_wind, ls_wind_stencil_extent, wind, m3_inv )

      USE constants_mod, ONLY: r_def, i_def
      USE field_mod, ONLY: field_type, field_proxy_type
      USE operator_mod, ONLY: operator_type, operator_proxy_type
      USE atl_w3h_advective_update_kernel_mod, ONLY: atl_w3h_advective_update_code
      USE mesh_mod, ONLY: mesh_type
      USE stencil_2D_dofmap_mod, ONLY: stencil_2D_dofmap_type, STENCIL_2D_CROSS

      TYPE(field_type), intent(in) :: advective_increment, tracer, ls_wind, wind
      TYPE(operator_type), intent(in) :: m3_inv
      INTEGER(KIND=i_def), intent(in) :: tracer_stencil_extent, ls_wind_stencil_extent
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:,:,:) :: m3_inv_local_stencil
      TYPE(operator_proxy_type) m3_inv_proxy
      REAL(KIND=r_def), pointer, dimension(:) :: wind_data
      REAL(KIND=r_def), pointer, dimension(:) :: ls_wind_data
      REAL(KIND=r_def), pointer, dimension(:) :: tracer_data
      REAL(KIND=r_def), pointer, dimension(:) :: advective_increment_data
      TYPE(field_proxy_type) advective_increment_proxy, tracer_proxy, ls_wind_proxy, wind_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_tracer(:,:), map_any_w2(:,:), map_w3(:,:)
      INTEGER(KIND=i_def) ndf_w3, undf_w3, ndf_adspc1_tracer, undf_adspc1_tracer, ndf_any_w2, undf_any_w2
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh
      INTEGER(KIND=i_def) ls_wind_max_branch_length
      INTEGER(KIND=i_def), pointer :: ls_wind_stencil_size(:,:)
      INTEGER(KIND=i_def), pointer :: ls_wind_stencil_dofmap(:,:,:,:)
      TYPE(stencil_2D_dofmap_type), pointer :: ls_wind_stencil_map
      INTEGER(KIND=i_def) tracer_max_branch_length
      INTEGER(KIND=i_def), pointer :: tracer_stencil_size(:,:)
      INTEGER(KIND=i_def), pointer :: tracer_stencil_dofmap(:,:,:,:)
      TYPE(stencil_2D_dofmap_type), pointer :: tracer_stencil_map
      !
      ! Initialise field and/or operator proxies
      !
      advective_increment_proxy = advective_increment%get_proxy()
      advective_increment_data => advective_increment_proxy%data
      tracer_proxy = tracer%get_proxy()
      tracer_data => tracer_proxy%data
      ls_wind_proxy = ls_wind%get_proxy()
      ls_wind_data => ls_wind_proxy%data
      wind_proxy = wind%get_proxy()
      wind_data => wind_proxy%data
      m3_inv_proxy = m3_inv%get_proxy()
      m3_inv_local_stencil => m3_inv_proxy%local_stencil
      !
      ! Initialise number of layers
      !
      nlayers = advective_increment_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => advective_increment_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Initialise stencil dofmaps
      !
      tracer_stencil_map => tracer_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,tracer_stencil_extent)
      tracer_max_branch_length = tracer_stencil_extent + 1_i_def
      tracer_stencil_dofmap => tracer_stencil_map%get_whole_dofmap()
      tracer_stencil_size => tracer_stencil_map%get_stencil_sizes()
      ls_wind_stencil_map => ls_wind_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,ls_wind_stencil_extent)
      ls_wind_max_branch_length = ls_wind_stencil_extent + 1_i_def
      ls_wind_stencil_dofmap => ls_wind_stencil_map%get_whole_dofmap()
      ls_wind_stencil_size => ls_wind_stencil_map%get_stencil_sizes()
      !
      ! Look-up dofmaps for each function space
      !
      map_w3 => advective_increment_proxy%vspace%get_whole_dofmap()
      map_adspc1_tracer => tracer_proxy%vspace%get_whole_dofmap()
      map_any_w2 => ls_wind_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for w3
      !
      ndf_w3 = advective_increment_proxy%vspace%get_ndf()
      undf_w3 = advective_increment_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc1_tracer
      !
      ndf_adspc1_tracer = tracer_proxy%vspace%get_ndf()
      undf_adspc1_tracer = tracer_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for any_w2
      !
      ndf_any_w2 = ls_wind_proxy%vspace%get_ndf()
      undf_any_w2 = ls_wind_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_halo_cell(1)
      !
      ! Call kernels and communication routines
      !
      IF (advective_increment_proxy%is_dirty(depth=1)) THEN
        CALL advective_increment_proxy%halo_exchange(depth=1)
      END IF
      !
      IF (tracer_proxy%is_dirty(depth=tracer_stencil_extent)) THEN
        CALL tracer_proxy%halo_exchange(depth=tracer_stencil_extent)
      END IF
      !
      IF (ls_wind_proxy%is_dirty(depth=ls_wind_stencil_extent)) THEN
        CALL ls_wind_proxy%halo_exchange(depth=ls_wind_stencil_extent)
      END IF
      !
      DO cell=loop0_start,loop0_stop
        !
        CALL atl_w3h_advective_update_code(cell, nlayers, advective_increment_data, tracer_data, tracer_stencil_size(:,cell), &
&tracer_max_branch_length, tracer_stencil_dofmap(:,:,:,cell), ls_wind_data, ls_wind_stencil_size(:,cell), &
&ls_wind_max_branch_length, ls_wind_stencil_dofmap(:,:,:,cell), wind_data, m3_inv_proxy%ncell_3d, m3_inv_local_stencil, ndf_w3, &
&undf_w3, map_w3(:,cell), ndf_adspc1_tracer, undf_adspc1_tracer, map_adspc1_tracer(:,cell), ndf_any_w2, undf_any_w2, &
&map_any_w2(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL advective_increment_proxy%set_dirty()
      CALL advective_increment_proxy%set_clean(1)
      CALL wind_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_atl_w3h_adv_upd_kernel_type

END MODULE invoke_atl_w3h_adv_upd_mod
