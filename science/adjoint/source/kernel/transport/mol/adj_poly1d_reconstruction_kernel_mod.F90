!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint of horizontal tracer values
!! through fitting a high order 1D upwind reconstruction.

! @todo Kernel does not work correctly when run in parallel; fix is being worked on

module adj_poly1d_reconstruction_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_READWRITE, GH_READ,       &
                              STENCIL, CROSS,              &
                              CELL_COLUMN,                 &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              ANY_DISCONTINUOUS_SPACE_2,   &
                              ANY_DISCONTINUOUS_SPACE_3
use constants_mod,     only : r_tran, i_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: adj_poly1d_reconstruction_kernel_type
type(arg_type) :: meta_args(6) = (/                                                            &
     arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! reconstruction
     arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),                 & ! tracer
     arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2, STENCIL(CROSS)), & ! dummy_ads2
     arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3),                 & ! coeff
     arg_type(gh_scalar, GH_INTEGER, GH_READ),                                                 & ! ndata
     arg_type(gh_scalar, GH_INTEGER, GH_READ)                                                  & ! order
     /)
integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: adj_poly1d_reconstruction_code
end type adj_poly1d_reconstruction_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_poly1d_reconstruction_code

contains
!> @brief Computes the horizontal polynomial interpolation of a tracer.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] reconstruction Reconstructed tracer field to compute
!> @param[in,out] tracer         Pointwise tracer field to reconstruct
!> @param[in]     stencil_size   Size of the stencil (number of cells)
!> @param[in]     stencil_map    Dofmaps for the stencil
!> @param[in]     coeff          Array of polynomial coefficients for interpolation
!> @param[in]     ndata          Number of data points per dof location
!> @param[in]     order          Desired order of polynomial reconstruction
!> @param[in]     ndf_md         Number of degrees of freedom per cell for reconstructed field
!> @param[in]     undf_md        Number of unique degrees of freedom for the
!!                               reconstructed field
!> @param[in]     map_md         Dofmap for the cell at the base of the column for reconstructed field
!> @param[in]     ndf_ws         Number of degrees of freedom per cell for the tracer
!> @param[in]     undf_ws        Number of unique degrees of freedom for the tracer field
!> @param[in]     map_ws         Dofmap for the cell at the base of the column for the tracer field
!> @param[in]     ndf_c          Number of degrees of freedom per cell for the coeff space
!> @param[in]     undf_c         Total number of degrees of freedom for the coeff space
!> @param[in]     map_c          Dofmap for the coeff space
subroutine adj_poly1d_reconstruction_code( nlayers, &
                                           reconstruction, &
                                           tracer, &
                                           dummy_ads2, &
                                           stencil_size, &
                                           stencil_map, &
                                           coeff, &
                                           ndata, &
                                           order, &
                                           ndf_md, &
                                           undf_md, &
                                           map_md, &
                                           ndf_ws, &
                                           undf_ws, &
                                           map_ws, &
                                           ndf_c, &
                                           undf_c, &
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
  integer(kind=i_def), intent(in)                                 :: order
  integer(kind=i_def), intent(in)                                 :: stencil_size
  real(kind=r_tran), dimension(undf_md), intent(inout)            :: reconstruction
  real(kind=r_tran), dimension(undf_ws), intent(inout)            :: tracer
  real(kind=r_tran), dimension(undf_ws), intent(in)               :: dummy_ads2
  real(kind=r_tran), dimension(undf_c), intent(in)                :: coeff
  integer(kind=i_def), dimension(ndf_ws,stencil_size), intent(in) :: stencil_map

  ! Internal variables
  integer(kind=i_def), parameter                      :: nfaces = 4
  integer(kind=i_def)                                 :: k
  integer(kind=i_def)                                 :: f
  integer(kind=i_def)                                 :: p
  integer(kind=i_def)                                 :: face
  integer(kind=i_def)                                 :: stencil
  integer(kind=i_def)                                 :: stencil_depth
  integer(kind=i_def)                                 :: depth
  integer(kind=i_def)                                 :: face_mod
  integer(kind=i_def)                                 :: ijp
  integer(kind=i_def)                                 :: df
  integer(kind=i_def)                                 :: nl
  integer(kind = i_def), dimension(order + 1, nfaces) :: map1d

  nl = ndf_ws + nlayers - 2
  stencil_depth = order / 2
  map1d(1,:) = 1
  do face = 1, nfaces, 1
    depth = 1
    face_mod = stencil_depth * MOD(face + 1, 2)
    do stencil = 2, stencil_depth + 1, 1
      map1d(stencil + depth - 1,face) = face_mod + stencil
      map1d(stencil + depth,face) = face_mod + order + stencil
      depth = depth + 1
    enddo
  enddo
  do f = nfaces, 1, -1
    df = f * nl + f - nl + map_md(1) - 1
    do p = order + 1, 1, -1
      ijp = f * order + f - order + p + map_c(1) - 2
      stencil = map1d(p,f)
      do k = nl, 0, -1
        tracer(k + stencil_map(1,stencil)) = tracer(k + stencil_map(1,stencil)) + coeff(ijp) * reconstruction(df + k)
      enddo
    enddo
    do k = nl, 0, -1
      reconstruction(df + k) = 0.0_r_tran
    enddo
  enddo

end subroutine adj_poly1d_reconstruction_code

end module adj_poly1d_reconstruction_kernel_mod
