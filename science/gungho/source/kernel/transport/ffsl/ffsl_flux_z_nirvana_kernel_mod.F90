!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes the vertical mass flux using Nirvana.
!> @details This kernel reconstructs a field using the reversible Nirvana scheme
!!          which is equivalent to fitting a quadratic to the cell such that the
!!          integral of the quadratic equals the integral of the field in the
!!          cell, and the quadratic matches the gradient of the field at
!!          cell edges. A limiter can be applied to ensure monotonicity.
!!          This kernel is designed to work in the vertical direction only and
!!          takes into account the vertical boundaries and grid spacing.
!!
!!          Note that this kernel only works when field is a W3 field at lowest
!!          order, since it is assumed that ndf_w3 = 1.

module ffsl_flux_z_nirvana_kernel_mod

use argument_mod,                   only : arg_type,              &
                                           GH_FIELD, GH_REAL,     &
                                           GH_READ, GH_WRITE,     &
                                           GH_SCALAR, GH_INTEGER, &
                                           GH_LOGICAL, CELL_COLUMN
use fs_continuity_mod,              only : W3, W2v
use constants_mod,                  only : r_tran, i_def, l_def, EPS_R_TRAN
use kernel_mod,                     only : kernel_type
use transport_enumerated_types_mod, only : vertical_monotone_relaxed, &
                                           vertical_monotone_strict,  &
                                           vertical_monotone_positive

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ffsl_flux_z_nirvana_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                  &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! flux
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! frac_wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! dz
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! detj
       arg_type(GH_SCALAR, GH_REAL,    GH_READ),       & ! dt
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),       & ! monotone
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)        & ! log_space
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: ffsl_flux_z_nirvana_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: ffsl_flux_z_nirvana_code

contains

!> @brief Compute the flux using the Nirvana reconstruction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] flux      The flux to be computed
!> @param[in]     frac_wind The fractional vertical wind
!> @param[in]     dep_dist  The vertical departure points
!> @param[in]     field     The field to construct the flux
!> @param[in]     dz        Vertical length of the W3 cell
!> @param[in]     detj      Volume of cells
!> @param[in]     dt        Time step
!> @param[in]     monotone  Monotonicity option to use
!> @param[in]     log_space Switch to use natural logarithmic space
!!                          for edge interpolation
!> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine ffsl_flux_z_nirvana_code( nlayers,   &
                                     flux,      &
                                     frac_wind, &
                                     dep_dist,  &
                                     field,     &
                                     dz,        &
                                     detj,      &
                                     dt,        &
                                     monotone,  &
                                     log_space, &
                                     ndf_w2v,   &
                                     undf_w2v,  &
                                     map_w2v,   &
                                     ndf_w3,    &
                                     undf_w3,   &
                                     map_w3 )

  use subgrid_horizontal_support_mod, only: bound_field
  use subgrid_vertical_support_mod,   only: vertical_ppm_recon,                  &
                                            vertical_ppm_mono_relax,             &
                                            vertical_ppm_mono_strict,            &
                                            vertical_ppm_positive,               &
                                            third_order_vertical_edge

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: ndf_w3
  real(kind=r_tran),   intent(inout) :: flux(undf_w2v)
  real(kind=r_tran),   intent(in)    :: field(undf_w3)
  real(kind=r_tran),   intent(in)    :: frac_wind(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dz(undf_w3)
  real(kind=r_tran),   intent(in)    :: detj(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  real(kind=r_tran),   intent(in)    :: dt
  integer(kind=i_def), intent(in)    :: monotone
  logical(kind=l_def), intent(in)    :: log_space

  ! Internal variables
  integer(kind=i_def) :: k, i, w2v_idx, w3_idx
  integer(kind=i_def) :: int_displacement, sign_displacement
  integer(kind=i_def) :: dep_cell_idx, sign_offset
  integer(kind=i_def) :: lowest_whole_cell, highest_whole_cell
  integer(kind=i_def) :: local_cell_idx

  real(kind=r_tran)   :: displacement, frac_dist
  real(kind=r_tran)   :: reconstruction, mass_whole_cells
  real(kind=r_tran)   :: edge_above(0:nlayers-1), edge_below(0:nlayers-1)
  real(kind=r_tran)   :: log_field(3)

  w3_idx = map_w3(1)
  w2v_idx = map_w2v(1)

  ! ========================================================================== !
  ! Reconstruct edge values
  ! ========================================================================== !

  ! Use Nirvana reconstruction for interpolated edge values

  if (log_space) then ! --------------------------------------------------------

    ! Bottom layer
    local_cell_idx = 1
    log_field(:) = LOG(MAX(ABS(field(w3_idx : w3_idx+2)), EPS_R_TRAN))
    call third_order_vertical_edge(log_field, dz(w3_idx : w3_idx+2),           &
                                   local_cell_idx, edge_above(0), edge_below(0))

    ! Centre layers
    local_cell_idx = 2
    do k = 1, nlayers - 2
      log_field(:) = LOG(MAX(ABS(field(w3_idx+k-1 : w3_idx+k+1)), EPS_R_TRAN))
      call third_order_vertical_edge(log_field, dz(w3_idx+k-1 : w3_idx+k+1),   &
                                     local_cell_idx, edge_above(k), edge_below(k))
    end do

    ! Top layer
    local_cell_idx = 3
    k = nlayers - 1
    log_field(:) = LOG(MAX(ABS(field(w3_idx+k-2 : w3_idx+k)), EPS_R_TRAN))
    call third_order_vertical_edge(log_field, dz(w3_idx+k-2 : w3_idx+k),       &
                                   local_cell_idx, edge_above(k), edge_below(k))

    ! Convert back from log space
    edge_above(:) = EXP(edge_above(:))
    edge_below(:) = EXP(edge_below(:))

  else ! Not log-space ---------------------------------------------------------

    ! Bottom layer
    local_cell_idx = 1
    call third_order_vertical_edge(field(w3_idx : w3_idx+2),                   &
                                   dz(w3_idx : w3_idx+2),                      &
                                   local_cell_idx, edge_above(0), edge_below(0))

    ! Centre layers
    local_cell_idx = 2
    do k = 1, nlayers - 2
      call third_order_vertical_edge(field(w3_idx+k-1 : w3_idx+k+1),           &
                                     dz(w3_idx+k-1 : w3_idx+k+1),              &
                                     local_cell_idx, edge_above(k), edge_below(k))
    end do

    ! Top layer
    local_cell_idx = 3
    k = nlayers - 1
    call third_order_vertical_edge(field(w3_idx+k-2 : w3_idx+k),               &
                                   dz(w3_idx+k-2 : w3_idx+k),                  &
                                   local_cell_idx, edge_above(k), edge_below(k))

  end if

  ! ========================================================================== !
  ! Apply monotonicity to edge reconstructions
  ! ========================================================================== !

  ! Loop again through edge reconstructions to limit them using neighbouring values

  select case ( monotone )
  case ( vertical_monotone_strict, vertical_monotone_relaxed )

    ! Bottom edge
    call bound_field( edge_below(0), field(w3_idx), field(w3_idx+1) )

    ! Centre edges
    do k = 1, nlayers-1
      call bound_field( edge_above(k-1), field(w3_idx+k-1), field(w3_idx+k) )
      call bound_field( edge_below(k), field(w3_idx+k-1), field(w3_idx+k) )
    end do

    ! Top edge
    call bound_field( edge_above(nlayers-1), field(w3_idx+nlayers-2),          &
                      field(w3_idx+nlayers-1) )

  case ( vertical_monotone_positive )

    ! Ensure all edges are positive
    do k = 0, nlayers-1
      edge_above(k) = max(edge_above(k),0.0_r_tran)
      edge_below(k) = max(edge_below(k),0.0_r_tran)
    end do

  end select

  ! ========================================================================== !
  ! Build fluxes
  ! ========================================================================== !

  ! Force bottom flux to be zero
  flux(w2v_idx) = 0.0_r_tran

  ! Loop through faces
  do k = 1, nlayers - 1

    ! Pull out departure point, and separate into integer / frac parts
    displacement = dep_dist(w2v_idx + k)
    int_displacement = INT(displacement, i_def)
    frac_dist = displacement - REAL(int_displacement, r_tran)
    sign_displacement = INT(SIGN(1.0_r_tran, displacement))

    ! Set an offset for the stencil index, based on dep point sign
    sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

    ! Determine departure cell
    dep_cell_idx = k - int_displacement + sign_offset - 1

    ! ======================================================================== !
    ! Integer sum
    mass_whole_cells = 0.0_r_tran
    lowest_whole_cell = MIN(dep_cell_idx + 1, k)
    highest_whole_cell = MAX(dep_cell_idx, k) - 1
    do i = lowest_whole_cell, highest_whole_cell
      mass_whole_cells = mass_whole_cells + field(w3_idx + i) * detj(w3_idx + i)
    end do

    ! ======================================================================== !
    ! Perform reversible PPM reconstruction for fractional part
    call vertical_ppm_recon(reconstruction, frac_dist,    &
                            field(w3_idx + dep_cell_idx), &  ! field in dep cell
                            edge_below(dep_cell_idx),     &  ! edge below dep cell
                            edge_above(dep_cell_idx))        ! edge above dep cell

    ! ======================================================================== !
    ! Apply monotonicity for fractional part
    select case ( monotone )
    case ( vertical_monotone_strict )
      call vertical_ppm_mono_strict(reconstruction,                            &
                                    field(w3_idx + dep_cell_idx),              &
                                    edge_below(dep_cell_idx),                  &
                                    edge_above(dep_cell_idx))

    case ( vertical_monotone_relaxed )
      call vertical_ppm_mono_relax(reconstruction, frac_dist,                  &
                                   field(w3_idx + dep_cell_idx),               &
                                   edge_below(dep_cell_idx),                   &
                                   edge_above(dep_cell_idx))

    case ( vertical_monotone_positive )
      call vertical_ppm_positive(reconstruction, frac_dist,                    &
                                 field(w3_idx + dep_cell_idx),                 &
                                 edge_below(dep_cell_idx),                     &
                                 edge_above(dep_cell_idx))

    end select

    ! ======================================================================== !
    ! Compute flux
    flux(w2v_idx + k) = (frac_wind(w2v_idx + k) * reconstruction               &
                         + sign(1.0_r_tran, displacement) * mass_whole_cells) / dt

  end do

  ! Force top flux to be zero
  flux(w2v_idx + nlayers) = 0.0_r_tran

end subroutine ffsl_flux_z_nirvana_code

end module ffsl_flux_z_nirvana_kernel_mod
