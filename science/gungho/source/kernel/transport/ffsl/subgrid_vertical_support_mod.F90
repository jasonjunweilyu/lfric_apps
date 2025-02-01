!------------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Routines for computing vertical Nirvana and PPM reconstructions
!!          of fields, for use in FFSL
!------------------------------------------------------------------------------
module subgrid_vertical_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN
use subgrid_horizontal_support_mod, only: bound_field

implicit none

private

! Edge reconstructions
public :: second_order_vertical_edge
public :: third_order_vertical_edge
public :: fourth_order_vertical_edge
! Field reconstructions
public :: vertical_ppm_recon
! Monotonic limiters
public :: fourth_order_vertical_mono
public :: fourth_order_vertical_quasi_mono
public :: vertical_ppm_mono_strict
public :: vertical_ppm_mono_relax
public :: vertical_ppm_positive

contains

  ! ========================================================================== !
  ! EDGE RECONSTRUCTIONS
  ! ========================================================================== !

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a second-order interpolation.
  !> @details Uses a second-order interpolation to find the vertical cell edge
  !!          values of the field. The vertical grid spacing is used to compute the
  !!          mass, and a quadratic is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge field value.
  !!
  !> @param[in]   field      Field values of two cells which have the ordering
  !!                         | 1 | 2 |
  !> @param[in]   dz         Height of each layer, with index the same as field
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 |
  !!                         with edges  0   1   2
  !> @param[out]  edge_value The interpolated field value at edge_to_do
  !----------------------------------------------------------------------------
  subroutine second_order_vertical_edge(field, dz, edge_to_do, edge_value)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(in)  :: field(2)
    real(kind=r_tran),   intent(in)  :: dz(2)
    integer(kind=i_def), intent(in)  :: edge_to_do
    real(kind=r_tran),   intent(out) :: edge_value

    ! Internal Variables
    real(kind=r_tran) :: z(0:2), edge_height
    real(kind=r_tran) :: cmass(0:2)

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    z(1) = z(0) + dz(1)
    z(2) = z(1) + dz(2)

    ! Get edge height to interpolate field
    edge_height = z(edge_to_do)

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    cmass(1) = cmass(0) + dz(1)*field(1)
    cmass(2) = cmass(1) + dz(2)*field(2)

    ! Calculate derivative of the quadratic at z = edge_height
    edge_value =   ( 2.0_r_tran*edge_height - z(2) ) / ( z(1) * ( z(1)-z(2) ) ) * cmass(1) &
                 + ( 2.0_r_tran*edge_height - z(1) ) / ( z(2) * ( z(2)-z(1) ) ) * cmass(2)

  end subroutine second_order_vertical_edge

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a third-order interpolation.
  !> @details Uses a third-order interpolation to find the vertical cell edge
  !!          values of the field. The vertical grid spacing is used to compute the
  !!          mass, and a high-order polynomial is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge value.
  !!          Both the top and bottom edges are computed for the cell, as this
  !!          edge is not continuous across cells.
  !!
  !> @param[in]   field      Field values of three cells which have the ordering
  !!                         | 1 | 2 | 3 |
  !> @param[in]   dz         Height of each layer, with index the same as field
  !> @param[in]   cell_to_do Tells routine which cell edges to do based on
  !!                         cells       | 1 | 2 | 3 |
  !!                         with edges  0   1   2   3
  !> @param[out]  edge_above The interpolated edge value located above the cell
  !> @param[out]  edge_below The interpolated edge value located below the cell
  !----------------------------------------------------------------------------
  subroutine third_order_vertical_edge(field, dz, cell_to_do, edge_above, edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: field(1:3)
    real(kind=r_tran),    intent(in)    :: dz(1:3)
    integer(kind=i_def),  intent(in)    :: cell_to_do
    real(kind=r_tran),    intent(out)   :: edge_above
    real(kind=r_tran),    intent(out)   :: edge_below

    real(kind=r_tran) :: z(0:3), dzs(1:3), dzsum, edge_height
    real(kind=r_tran) :: dmass(1:3)
    real(kind=r_tran) :: cmass(0:3)
    real(kind=r_tran) :: poly_mass(1:3)
    real(kind=r_tran) :: dl_dz(1:3)

    integer(kind=i_def) :: i

    ! Get scaling value
    dzsum = sum(dz)

    ! Get scaled dz
    dzs = dz / dzsum

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    do i = 1, 3
      z(i) = z(i-1) + dzs(i)
    end do

    ! Get mass scaled by height
    dmass = field * dzs

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    do i = 1, 3
      cmass(i) = cmass(i-1) + dmass(i)
    end do

    ! Get cumulative mass divided by denominator of polynomial
    poly_mass(1) = cmass(1) / ( (z(1))*(z(1)-z(2))*(z(1)-z(3)) )
    poly_mass(2) = cmass(2) / ( (z(2))*(z(2)-z(1))*(z(2)-z(3)) )
    poly_mass(3) = cmass(3) / ( (z(3))*(z(3)-z(1))*(z(3)-z(2)) )

    ! Get edge height to interpolate field to for edge below the cell
    edge_height = z(cell_to_do-1)

    ! Calculate derivative of numerator of polynomial at edge height
    dl_dz(1) = 3.0_r_tran*edge_height**2 &
               - 2.0_r_tran*(z(2)+z(3))*edge_height + z(2)*z(3)
    dl_dz(2) = 3.0_r_tran*edge_height**2 &
               - 2.0_r_tran*(z(1)+z(3))*edge_height + z(1)*z(3)
    dl_dz(3) = 3.0_r_tran*edge_height**2 &
               - 2.0_r_tran*(z(1)+z(2))*edge_height + z(1)*z(2)

    ! Calculate value of edge below cell
    edge_below = sum( poly_mass * dl_dz )

    ! Get edge height to interpolate field to for edge above the cell
    edge_height = z(cell_to_do)

    ! Calculate derivative of numerator of polynomial at edge height
    dl_dz(1) = 3.0_r_tran*edge_height**2 &
               - 2.0_r_tran*(z(2)+z(3))*edge_height + z(2)*z(3)
    dl_dz(2) = 3.0_r_tran*edge_height**2 &
               - 2.0_r_tran*(z(1)+z(3))*edge_height + z(1)*z(3)
    dl_dz(3) = 3.0_r_tran*edge_height**2 &
               - 2.0_r_tran*(z(1)+z(2))*edge_height + z(1)*z(2)

    ! Calculate value of edge above cell
    edge_above = sum( poly_mass * dl_dz )

  end subroutine third_order_vertical_edge

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a fourth-order interpolation.
  !> @details Uses a fourth-order interpolation to find the vertical cell edge
  !!          values of the field. The vertical grid spacing is used to compute the
  !!          mass, and a high-order polynomial is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge value.
  !!
  !> @param[in]   field      Field values of four cells which have the ordering
  !!                         | 1 | 2 | 3 | 4 |
  !> @param[in]   dz         Height of each layer, with index the same as field
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 | 3 | 4 |
  !!                         with edges  0   1   2   3   4
  !> @param[out]  edge_below Interpolated field value at edge_to_do
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_edge(field, dz, edge_to_do, edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: field(1:4)
    real(kind=r_tran),    intent(in)    :: dz(1:4)
    integer(kind=i_def),  intent(in)    :: edge_to_do
    real(kind=r_tran),    intent(out)   :: edge_below

    real(kind=r_tran) :: z(0:4), dzs(1:4), dzsum, edge_height
    real(kind=r_tran) :: dmass(1:4)
    real(kind=r_tran) :: cmass(0:4)
    real(kind=r_tran) :: poly_mass(1:4)
    real(kind=r_tran) :: dl_dz(1:4)

    integer(kind=i_def) :: i

    ! Get scaling value
    dzsum = sum(dz)

    ! Get scaled dz
    dzs = dz / dzsum

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    do i = 1, 4
      z(i) = z(i-1) + dzs(i)
    end do

    ! Get edge height to interpolate field to
    edge_height = z(edge_to_do)

    ! Get mass scaled by height
    dmass = field * dzs

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    do i = 1, 4
      cmass(i) = cmass(i-1) + dmass(i)
    end do

    ! Get cumulative mass divided by denominator of polynomial
    poly_mass(1) = cmass(1)/((z(1))*(z(1)-z(2))*(z(1)-z(3))*(z(1)-z(4)))
    poly_mass(2) = cmass(2)/((z(2))*(z(2)-z(1))*(z(2)-z(3))*(z(2)-z(4)))
    poly_mass(3) = cmass(3)/((z(3))*(z(3)-z(1))*(z(3)-z(2))*(z(3)-z(4)))
    poly_mass(4) = cmass(4)/((z(4))*(z(4)-z(1))*(z(4)-z(2))*(z(4)-z(3)))

    ! Calculate derivative of numerator of polynomial at edge height
    dl_dz    = 4.0_r_tran*edge_height**3
    dl_dz(1) = dl_dz(1) - 3.0_r_tran*(z(2)+z(3)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(3)*z(4) + z(2)*z(3) + z(2)*z(4))*edge_height - z(2)*z(3)*z(4)
    dl_dz(2) = dl_dz(2) - 3.0_r_tran*(z(1)+z(3)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(3)*z(4) + z(1)*z(3) + z(1)*z(4))*edge_height - z(1)*z(3)*z(4)
    dl_dz(3) = dl_dz(3) - 3.0_r_tran*(z(1)+z(2)+z(4))*edge_height**2 &
               + 2.0_r_tran*(z(2)*z(4) + z(1)*z(2) + z(1)*z(4))*edge_height - z(1)*z(2)*z(4)
    dl_dz(4) = dl_dz(4) - 3.0_r_tran*(z(1)+z(2)+z(3))*edge_height**2 &
               + 2.0_r_tran*(z(2)*z(3) + z(1)*z(2) + z(1)*z(3))*edge_height - z(1)*z(2)*z(3)

    ! Calculate value of edge below layer k
    edge_below = sum( poly_mass * dl_dz )

  end subroutine fourth_order_vertical_edge

  ! ========================================================================== !
  ! FIELD RECONSTRUCTIONS
  ! ========================================================================== !

  !----------------------------------------------------------------------------
  !> @brief  Returns the vertical PPM reconstruction (also used for the reversible
  !!         Nirvana reconstruction). This can be used to compute the flux as:
  !!         flux = u * reconstruction
  !!         The reconstruction is a third-order reconstruction at the cell's
  !!         vertical edges, and is based on a quadratic subgrid reconstruction that
  !!         uses the field interpolated to cell edges.
  !!
  !! @param[out]  recon       The PPM reconstruction
  !! @param[in]   dep         The fractional departure distance for the reconstruction point.
  !!                          For dep>0 the recon lives above the cell
  !!                          For dep<0 the recon lives below the cell
  !! @param[in]   field       Field values in the cell
  !! @param[in]   edge_below  Estimate of edge value below the cell
  !! @param[in]   edge_above  Estimate of edge value above the cell
  !----------------------------------------------------------------------------
  subroutine vertical_ppm_recon(recon,      &
                                dep,        &
                                field,      &
                                edge_below, &
                                edge_above)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(out) :: recon
    real(kind=r_tran),   intent(in)  :: dep
    real(kind=r_tran),   intent(in)  :: field
    real(kind=r_tran),   intent(in)  :: edge_below
    real(kind=r_tran),   intent(in)  :: edge_above

    ! Reconstruction weights
    real(kind=r_tran) :: cm, cc, cp

    ! Compute reconstruction weights
    if (dep >= 0.0_r_tran) then
      cp = 1.0_r_tran - 2.0_r_tran*dep + dep**2
      cc = 3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = - dep + dep**2
    else
      cp = dep + dep**2
      cc = -3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = 1.0_r_tran + 2.0_r_tran*dep + dep**2
    end if

    ! Apply weights to field and field edge values
    recon = cm*edge_below + cc*field + cp*edge_above

  end subroutine vertical_ppm_recon

  ! ========================================================================== !
  ! MONOTONIC LIMITERS
  ! ========================================================================== !

  !----------------------------------------------------------------------------
  !> @brief Applies monotonicity to a 4th-order edge reconstruction.
  !> @param[in]    field      Field values of four cells which have the ordering
  !!                          | 1 | 2 | 3 | 4 |
  !> @param[in]    edge_to_do Tells routine which edge to do based on
  !!                          cells       | 1 | 2 | 3 | 4 |
  !!                          with edges  0   1   2   3   4
  !> @param[inout] edge_below The edge value located at edge_to_do
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_mono(field,      &
                                        edge_to_do, &
                                        edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: field(1:4)
    integer(kind=i_def),  intent(in)    :: edge_to_do
    real(kind=r_tran),    intent(inout) :: edge_below

    real(kind=r_tran) :: t1

    ! Strict Monotonicity
    if ( edge_to_do > 0_i_def .AND. edge_to_do < 4_i_def) then
      t1 = ( edge_below - field(edge_to_do) )*( field(edge_to_do+1) - edge_below )
      if ( t1 < 0.0_r_tran ) then
        call bound_field(edge_below, field(edge_to_do), field(edge_to_do+1))
      end if
    else if ( edge_to_do == 0_i_def ) then
      call bound_field(edge_below, field(1), field(2))
    else if ( edge_to_do == 4_i_def ) then
      call bound_field(edge_below, field(3), field(4))
    end if

  end subroutine fourth_order_vertical_mono

  !----------------------------------------------------------------------------
  !> @brief Applies quasi-monotonic positivity to a 4th-order edge reconstruction.
  !!        This requires two cells either side of the edge.
  !> @param[in]    field      Field values of four cells which have the ordering
  !!                          | 1 | 2 | 3 | 4 |
  !! @param[in]    dep        The fractional departure distance at the edge
  !> @param[in]    min_val    Minimum value to enforce edge value to be
  !> @param[inout] edge_below The edge value located at edge 2 below
  !!                          cells       | 1 | 2 | 3 | 4 |
  !!                          with edges  0   1   2   3   4
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_quasi_mono(field,      &
                                              dep,        &
                                              min_val,    &
                                              edge_below)

    implicit none

    real(kind=r_tran),    intent(in)    :: field(1:4)
    real(kind=r_tran),    intent(in)    :: dep
    real(kind=r_tran),    intent(in)    :: min_val
    real(kind=r_tran),    intent(inout) :: edge_below

    real(kind=r_tran)   :: t1
    integer(kind=i_def) :: sign_dep, sign_cor, cell, cell_up, cell_dw

    ! Quasi-monotonic positive edges

    ! Get sign of departure distance
    sign_dep = INT(SIGN(1.0_r_tran, dep))
    sign_cor = INT((1.0_r_tran - SIGN(1.0_r_tran, dep))/2.0_r_tran)
    cell    = 2+sign_cor
    cell_up = 2+sign_cor+sign_dep
    cell_dw = 2+sign_cor-sign_dep
    ! Look at sign of successive gradients at edge in upwind direction
    t1 = ( field(cell_up) - field(cell) )*( field(cell) - field(cell_dw) )
    if ( t1 < 0.0_r_tran ) then
      call bound_field(edge_below, field(2), field(3))
    end if
    ! Apply positivity
    edge_below = max(edge_below, min_val)

  end subroutine fourth_order_vertical_quasi_mono

  !----------------------------------------------------------------------------
  !> @brief  Applies a strict limiter to a PPM reconstruction.
  !!
  !! @param[inout] recon       The PPM reconstruction
  !! @param[in]    field       Field values in the cell
  !! @param[in]    edge_below  Estimate of edge value below the cell
  !! @param[in]    edge_above  Estimate of edge value above the cell
  !----------------------------------------------------------------------------
  subroutine vertical_ppm_mono_strict(recon,      &
                                      field,      &
                                      edge_below, &
                                      edge_above)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: field
    real(kind=r_tran),   intent(in)    :: edge_below
    real(kind=r_tran),   intent(in)    :: edge_above

    ! Monotonicity variable
    real(kind=r_tran) :: t1

    ! Strict monotonicity
    t1 = (2.0_r_tran*edge_below + edge_above - 3.0_r_tran*field) &
         / (3.0_r_tran*edge_below + 3.0_r_tran*edge_above - 6.0_r_tran*field + EPS_R_TRAN)
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
      ! If subgrid reconstruction has extrema in the cell then revert to constant reconstruction
      recon = field
    end if

  end subroutine vertical_ppm_mono_strict

  !----------------------------------------------------------------------------
  !> @brief  Applies a relaxed limiter to a vertical PPM reconstruction.
  !!
  !! @param[inout] recon       The PPM reconstruction
  !! @param[in]    dep         The fractional departure distance for the reconstruction point.
  !!                           For dep>0 the recon lives above the cell
  !!                           For dep<0 the recon lives below the cell
  !! @param[in]    field       Field values in the cell
  !! @param[in]    edge_below  Estimate of edge value below the cell
  !! @param[in]    edge_above  Estimate of edge value above the cell
  !----------------------------------------------------------------------------
  subroutine vertical_ppm_mono_relax(recon,      &
                                     dep,        &
                                     field,      &
                                     edge_below, &
                                     edge_above)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: dep
    real(kind=r_tran),   intent(in)    :: field
    real(kind=r_tran),   intent(in)    :: edge_below
    real(kind=r_tran),   intent(in)    :: edge_above

    ! Internal variables
    real(kind=r_tran) :: cm, cc, cp
    real(kind=r_tran) :: t1, t2, t3

    ! Relaxed monotonicity
    t1 = (2.0_r_tran*edge_below + edge_above - 3.0_r_tran*field) &
         / (3.0_r_tran*edge_below + 3.0_r_tran*edge_above - 6.0_r_tran*field + EPS_R_TRAN)
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
      ! If subgrid reconstruction has extrema in the cell then check smoothness of field
      t2 = (edge_above - field)*(field - edge_below)
      t3 = abs(field - edge_below) - abs(edge_above - field)
      if ( t2 < 0.0_r_tran ) then
        ! Revert to constant reconstruction
        recon = field
      else if ( t3 < 0.0_r_tran .and. dep >= 0.0_r_tran ) then
        ! Ensure subgrid reconstruction is bounded by edge values
        cp = 0.0_r_tran
        cc = 3.0_r_tran - 3.0_r_tran*dep + dep**2
        cm = -2.0_r_tran + 3.0_r_tran*dep - dep**2
        recon = cm*edge_below + cc*field + cp*edge_above
      else if ( t3 < 0.0_r_tran ) then
        cp = 0.0_r_tran
        cc = dep**2
        cm = 1.0_r_tran - dep**2
        recon = cm*edge_below + cc*field + cp*edge_above
      else if (dep >= 0.0_r_tran) then
        cp = 1.0_r_tran - dep**2
        cc = dep**2
        cm = 0.0_r_tran
        recon = cm*edge_below + cc*field + cp*edge_above
      else
        cp = -2.0_r_tran - 3.0_r_tran*dep - dep**2
        cc = 3.0_r_tran + 3.0_r_tran*dep + dep**2
        cm = 0.0_r_tran
        recon = cm*edge_below + cc*field + cp*edge_above
      end if
    end if

  end subroutine vertical_ppm_mono_relax

  !----------------------------------------------------------------------------
  !> @brief  Applies a positive limiter to a vertical PPM reconstruction.
  !!
  !! @param[inout] recon       The PPM reconstruction
  !! @param[in]    dep         The fractional departure distance for the reconstruction point.
  !!                           For dep>0 the recon lives above the cell
  !!                           For dep<0 the recon lives below the cell
  !! @param[in]    field       Field values in the cell
  !! @param[in]    edge_below  Estimate of edge value below the cell
  !! @param[in]    edge_above  Estimate of edge value above the cell
  !----------------------------------------------------------------------------
  subroutine vertical_ppm_positive(recon,      &
                                   dep,        &
                                   field,      &
                                   edge_below, &
                                   edge_above)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(inout) :: recon
    real(kind=r_tran),   intent(in)    :: dep
    real(kind=r_tran),   intent(in)    :: field
    real(kind=r_tran),   intent(in)    :: edge_below
    real(kind=r_tran),   intent(in)    :: edge_above

    ! Internal variables
    real(kind=r_tran) :: t1, t2, aa, bb

    ! Positive definite limiter
    aa = -4.0_r_tran*edge_below - 2.0_r_tran*edge_above + 6.0_r_tran*field
    bb = 3.0_r_tran*edge_below + 3.0_r_tran*edge_above - 6.0_r_tran*field
    ! Find stationary point of parabolic subgrid reconstruction
    t1 = -0.5_r_tran*aa/(bb+EPS_R_TRAN)
    if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
      ! If stationary point lies within the grid cell and makes the subgrid reconstruction
      ! negative, we revert to constant reconstruction
      t2 = edge_below + aa * t1 + bb  * t1 * t1
      if ( t2 < 0.0_r_tran ) then
        recon = field
      end if
    end if

  end subroutine vertical_ppm_positive

end module subgrid_vertical_support_mod
