!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the adjoint miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module adjoint_test_driver_mod

  use base_mesh_config_mod,                    only : prime_mesh_name
  use field_mod,                               only : field_type
  use log_mod,                                 only : log_event, LOG_LEVEL_INFO
  use mesh_mod,                                only : mesh_type
  use mesh_collection_mod,                     only : mesh_collection
  use model_clock_mod,                         only : model_clock_type
  use sci_geometric_constants_mod,             only : get_coordinates, &
                                                      get_panel_id
  use setup_test_tl_transport_runtime_alg_mod, only : setup_test_tl_transport_runtime

  implicit none

  private

  public run

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief   Runs adjoint tests.
  !> @details Runs algorithm layer adjoint tests.
  subroutine run( model_clock )

    ! PSyAD generated tests
    use gen_adj_kernel_tests_mod,           only : run_gen_adj_kernel_tests

    ! Handwritten tests
    ! transport/common
    use adjt_end_transport_step_alg_mod,    only : adjt_build_up_flux_alg
    use atlt_end_transport_step_alg_mod,    only : atlt_end_adv_step_alg, &
                                                   atlt_end_con_step_alg
    ! transport/control
    use atlt_transport_field_alg_mod,       only : atlt_transport_field_alg
    ! transport/mol
    use atlt_poly_adv_update_alg_mod,       only : atlt_poly_adv_update_alg
    use adjt_poly_adv_update_alg_mod,       only : adjt_poly_adv_update_alg
    use adjt_poly1d_reconstruction_alg_mod, only : adjt_poly1d_reconstruction_alg
    use adjt_poly2d_reconstruction_alg_mod, only : adjt_poly2d_reconstruction_alg
    use atlt_poly1d_vert_w3_recon_alg_mod,  only : atlt_poly1d_vert_w3_recon_alg
    use adjt_w3h_advective_update_alg_mod,  only : adjt_w3h_advective_update_alg
    use atlt_w3h_advective_update_alg_mod,  only : atlt_w3h_advective_update_alg
    use atlt_wt_advective_update_alg_mod,   only : atlt_hori_wt_update_alg, &
                                                   atlt_vert_wt_update_alg, &
                                                   atlt_wt_advective_update_alg
    use adjt_wt_advective_update_alg_mod,   only : adjt_hori_wt_update_alg
    use atlt_reconstruct_w3_field_alg_mod,  only : atlt_vert_w3_reconstruct_alg, &
                                                   atlt_reconstruct_w3_field_alg
    use adjt_reconstruct_w3_field_alg_mod,  only : adjt_hori_w3_reconstruct_alg
    use atlt_advective_and_flux_alg_mod,    only : atlt_advective_and_flux_alg
    use atlt_mol_advective_alg_mod,         only : atlt_mol_advective_alg
    use atlt_mol_conservative_alg_mod,      only : atlt_mol_conservative_alg
    ! misc
    use adjt_interpolation_alg_mod,         only : adjt_interp_w3wth_to_w2_alg, &
                                                   adjt_interp_w2_to_w3wth_alg
    use adjt_convert_cart2sphere_vector_alg_mod, &
                                            only : adjt_convert_cart2sphere_vector_alg
    use adjt_sci_convert_hdiv_field_alg_mod, &
                                            only : adjt_sci_convert_hdiv_field_alg

    implicit none

    ! Arguments
    type(model_clock_type), intent(in) :: model_clock

    ! Internal variables
    type(field_type),          pointer :: chi(:)
    type(field_type),          pointer :: panel_id
    type(mesh_type),           pointer :: mesh

    mesh => mesh_collection%get_mesh( prime_mesh_name )
    chi => get_coordinates( mesh%get_id() )
    panel_id => get_panel_id( mesh%get_id() )

    call setup_test_tl_transport_runtime( mesh, model_clock )

    call log_event( "TESTING generated adjoint kernels", LOG_LEVEL_INFO )
    call run_gen_adj_kernel_tests( mesh, chi, panel_id )

    call log_event( "TESTING adjoint kernels", LOG_LEVEL_INFO )
    ! transport/mol
    call atlt_poly_adv_update_alg( mesh )
    call adjt_poly_adv_update_alg( mesh )
    call adjt_poly1d_reconstruction_alg( mesh )
    call adjt_poly2d_reconstruction_alg( mesh )
    call atlt_poly1d_vert_w3_recon_alg( mesh )
    call adjt_w3h_advective_update_alg( mesh )
    call atlt_w3h_advective_update_alg( mesh )

    call log_event( "TESTING adjoint algorithms", LOG_LEVEL_INFO )
    ! misc
    call adjt_convert_cart2sphere_vector_alg( mesh )
    call adjt_sci_convert_hdiv_field_alg( mesh, chi, panel_id )
    call adjt_interp_w3wth_to_w2_alg( mesh )
    call adjt_interp_w2_to_w3wth_alg( mesh )
    ! transport/common
    call adjt_build_up_flux_alg( mesh )
    call atlt_end_adv_step_alg( mesh )
    call atlt_end_con_step_alg( mesh )
    ! transport/mol
    call atlt_hori_wt_update_alg( mesh )
    call atlt_vert_wt_update_alg( mesh )
    call adjt_hori_wt_update_alg( mesh )
    call atlt_wt_advective_update_alg( mesh, model_clock )
    call atlt_vert_w3_reconstruct_alg( mesh )
    call adjt_hori_w3_reconstruct_alg( mesh )
    call atlt_reconstruct_w3_field_alg( mesh )
    call atlt_advective_and_flux_alg( mesh, model_clock )
    call atlt_mol_advective_alg( mesh )
    call atlt_mol_conservative_alg( mesh )
    ! transport/control
    call atlt_transport_field_alg( mesh, model_clock )

    call log_event( "TESTING COMPLETE", LOG_LEVEL_INFO )

  end subroutine run

end module adjoint_test_driver_mod
