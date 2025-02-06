!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint of horizontal tracer values
!!        through fitting a high order 2D upwind reconstruction.

! @todo Kernel does not work correctly when run in parallel; fix is being worked on

module adj_poly2d_reconstruction_kernel_mod

use argument_mod,      only : arg_type, func_type,       &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_READWRITE, GH_READ,     &
                              STENCIL, REGION,           &
                              CELL_COLUMN,               &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3
use constants_mod,     only : r_tran, i_def, l_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: adj_poly2d_reconstruction_kernel_type
  type(arg_type) :: meta_args(6) = (/ &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),                  & ! reconstruction
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),                  & ! tracer
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2, STENCIL(REGION)), & ! dummy_ads2
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3),                  & ! coeff
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                                  & ! ndata
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                                                   & ! stencil_size
       /)
  integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_poly2d_reconstruction_code
end type adj_poly2d_reconstruction_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_poly2d_reconstruction_code

contains

!> @brief Computes the adjoint of horizontal polynomial interpolation of a tracer.
!> @param[in]     nlayers          Number of layers
!> @param[in,out] reconstruction   Reconstructed tracer field to compute
!> @param[in,out] tracer           Pointwise tracer field to reconstruct
!> @param[in]     dummy_adcsp2     Dummy tracer field for stencil
!> @param[in]     cells_in_stencil Number of cells needed to compute the polynomial
!!                                 fit (may be less than stencil_size)
!> @param[in]     stencil_map      Dofmaps for the stencil
!> @param[in]     coeff            Array of polynomial coefficients for interpolation
!> @param[in]     ndata            Number of data points per dof location
!> @param[in]     stencil_size     Size of the stencil (number of cells)
!> @param[in]     ndf_md           Number of degrees of freedom per cell
!> @param[in]     undf_md          Number of unique degrees of freedom for the
!!                                 reconstructed field
!> @param[in]     map_md           Dofmap for the cell at the base of the column
!> @param[in]     ndf_ws           Number of degrees of freedom per cell
!> @param[in]     undf_ws          Number of unique degrees of freedom for the tracer field
!> @param[in]     map_ws           Dofmap for the cell at the base of the column for the tracer field
!> @param[in]     ndf_c            Number of degrees of freedom per cell for the coeff space
!> @param[in]     undf_c           Total number of degrees of freedom for the coeff space
!> @param[in]     map_c            Dofmap for the coeff space
subroutine adj_poly2d_reconstruction_code( nlayers,          &
                                           reconstruction,   &
                                           tracer,           &
                                           dummy_ads2,       &
                                           cells_in_stencil, &
                                           stencil_map,      &
                                           coeff,            &
                                           ndata,            &
                                           stencil_size,     &
                                           ndf_md,           &
                                           undf_md,          &
                                           map_md,           &
                                           ndf_ws,           &
                                           undf_ws,          &
                                           map_ws,           &
                                           ndf_c,            &
                                           undf_c,           &
                                           map_c )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                                 :: nlayers
  integer(kind=i_def), intent(in)                                 :: ndf_ws
  integer(kind=i_def), intent(in)                                 :: undf_ws
  integer(kind=i_def), dimension(ndf_ws), intent(in)              :: map_ws
  integer(kind=i_def), intent(in)                                 :: ndf_md
  integer(kind=i_def), intent(in)                                 :: undf_md
  integer(kind=i_def), dimension(ndf_md), intent(in)              :: map_md
  integer(kind=i_def), intent(in)                                 :: ndf_c
  integer(kind=i_def), intent(in)                                 :: undf_c
  integer(kind=i_def), dimension(ndf_c), intent(in)               :: map_c
  integer(kind=i_def), intent(in)                                 :: ndata
  integer(kind=i_def), intent(in)                                 :: cells_in_stencil
  integer(kind=i_def), intent(in)                                 :: stencil_size
  real(kind=r_tran), dimension(undf_md), intent(inout)            :: reconstruction
  real(kind=r_tran), dimension(undf_ws), intent(inout)            :: tracer
  real(kind=r_tran), dimension(undf_ws), intent(in)               :: dummy_ads2
  real(kind=r_tran), dimension(undf_c), intent(in)                :: coeff
  integer(kind=i_def), dimension(ndf_ws,stencil_size), intent(in) :: stencil_map

  ! Internal variables
  integer(kind=i_def), parameter :: nfaces = 4
  integer(kind=i_def)            :: k
  integer(kind=i_def)            :: f
  integer(kind=i_def)            :: p
  integer(kind=i_def)            :: id_r
  integer(kind=i_def)            :: id_t
  integer(kind=i_def)            :: id_c
  integer(kind=i_def)            :: nl

  nl = ndf_ws + nlayers - 2
  do f = nfaces, 1, -1
    id_r = f * nl + f - nl + map_md(1) - 1
    do p = cells_in_stencil, 1, -1
      id_c = f * stencil_size + p - stencil_size + map_c(1) - 1
      id_t = stencil_map(1,p)
      do k = nl, 0, -1
        tracer(id_t + k) = tracer(id_t + k) + coeff(id_c) * reconstruction(id_r + k)
      enddo
    enddo
    do k = nl, 0, -1
      reconstruction(id_r + k) = 0.0_r_tran
    enddo
  enddo

end subroutine adj_poly2d_reconstruction_code

end module adj_poly2d_reconstruction_kernel_mod
