!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes adjoint of horizontal advective update
!!        through fitting a high order upwind reconstruction.

! @todo Kernel does not work correctly when run in parallel; fix is being worked on

module adj_w3h_advective_update_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_REAL,         &
                              GH_OPERATOR,               &
                              GH_WRITE, GH_READ,     &
                              STENCIL, CROSS2D,          &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_W2, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def, l_def, r_tran
use fs_continuity_mod, only : W3
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: adj_w3h_advective_update_kernel_type
  type(arg_type) :: meta_args(5) = (/                                                             &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,      W3),                                          &
       arg_type(GH_FIELD,    GH_REAL, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),                   &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,      ANY_W2, STENCIL(CROSS2D)),                    &
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,      W3, W3)                                       &
       /)
  integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: adj_w3h_advective_update_code
end type adj_w3h_advective_update_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: adj_w3h_advective_update_code

contains
!> @brief Computes the adjoint horizontal advective update for a tracer in W3.
!> @param[in]     cell                Horizontal cell index
!> @param[in]     nlayers             Number of layers
!> @param[in]     advective_increment Advective update field to compute
!> @param[in,out] tracer              Pointwise tracer field to advect stored on cell faces
!> @param[in]     dummy_adcsp1        Read only tracer field for stencil
!> @param[in]     smap_md_size        Sizes of the stencil map in each direction
!> @param[in]     smap_md_max         Maximum size of the stencil map
!> @param[in]     smap_md             Stencil map for the tracer space
!> @param[in]     wind                Wind field
!> @param[in]     smap_w2_size        Sizes of the stencil map in each direction
!> @param[in]     smap_w2_max         Maximum size of the stencil map
!> @param[in]     smap_w2             Stencil map for the wind space
!> @param[in]     ncell_3d            Total number of cells
!> @param[in]     m3_inv              Inverse mass matrix for W3 space
!> @param[in]     ndf_w3              Number of degrees of freedom per cell
!> @param[in]     undf_w3             Number of unique degrees of freedom for the
!!                                    advective_update field
!> @param[in]     map_w3              Dofmap for the cell at the base of the column
!> @param[in]     ndf_wd              Number of degrees of freedom per cell
!> @param[in]     undf_wd             Number of unique degrees of freedom for the
!!                                    tracer field
!> @param[in]     map_wd              Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2              Number of degrees of freedom per cell for the wind fields
!> @param[in]     undf_w2             Number of unique degrees of freedom for the wind fields
!> @param[in]     map_w2              Dofmap for the cell at the base of the column for the wind fields
subroutine adj_w3h_advective_update_code( cell,                &
                                          nlayers,             &
                                          advective_increment, &
                                          tracer,              &
                                          dummy_adcsp1,        &
                                          smap_md_size,        &
                                          smap_md_max,         &
                                          smap_md,             &
                                          wind,                &
                                          smap_w2_size,        &
                                          smap_w2_max,         &
                                          smap_w2,             &
                                          ncell_3d,            &
                                          m3_inv,              &
                                          ndf_w3,              &
                                          undf_w3,             &
                                          map_w3,              &
                                          ndf_md,              &
                                          undf_md,             &
                                          map_md,              &
                                          ndf_w2,              &
                                          undf_w2,             &
                                          map_w2 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                                  :: nlayers
  integer(kind=i_def), intent(in)                                  :: cell
  integer(kind=i_def), intent(in)                                  :: ncell_3d
  integer(kind=i_def), intent(in)                                  :: ndf_w3
  integer(kind=i_def), intent(in)                                  :: ndf_w2
  integer(kind=i_def), intent(in)                                  :: ndf_md
  integer(kind=i_def), intent(in)                                  :: undf_w3
  integer(kind=i_def), intent(in)                                  :: undf_w2
  integer(kind=i_def), intent(in)                                  :: undf_md
  integer(kind=i_def), dimension(ndf_w3), intent(in)               :: map_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in)               :: map_w2
  integer(kind=i_def), dimension(ndf_md), intent(in)               :: map_md
  integer(kind=i_def), intent(in)                                  :: smap_md_max
  integer(kind=i_def), dimension(4), intent(in)                    :: smap_md_size
  integer(kind=i_def), dimension(ndf_md,smap_md_max,4), intent(in) :: smap_md
  integer(kind=i_def), intent(in)                                  :: smap_w2_max
  integer(kind=i_def), dimension(4), intent(in)                    :: smap_w2_size
  integer(kind=i_def), dimension(ndf_w2,smap_w2_max,4), intent(in) :: smap_w2
  real(kind=r_tran), dimension(undf_w3), intent(in)                :: advective_increment
  real(kind=r_tran), dimension(undf_w2), intent(in)                :: wind
  real(kind=r_tran), dimension(undf_md), intent(inout)             :: tracer
  real(kind=r_tran), dimension(undf_md), intent(in)                :: dummy_adcsp1
  real(kind=r_def), dimension(ncell_3d,ndf_w3,ndf_w3), intent(in)  :: m3_inv

  ! Internal variables
  integer(kind=i_def), parameter         :: nfaces = 4
  integer(kind=i_def)                    :: k
  integer(kind=i_def)                    :: ik
  integer(kind=i_def)                    :: face
  integer(kind=i_def)                    :: df
  integer(kind=i_def)                    :: df1
  integer(kind=i_def)                    :: df2
  integer(kind=r_tran)                   :: i_N
  integer(kind=r_tran)                   :: i_S
  integer(kind=r_tran)                   :: i_E
  integer(kind=r_tran)                   :: i_W
  real(kind=r_tran)                      :: u
  real(kind=r_tran)                      :: v
  real(kind=r_tran)                      :: dtdx
  real(kind=r_tran)                      :: dtdy
  real(kind=r_tran)                      :: t_E
  real(kind=r_tran)                      :: t_W
  real(kind=r_tran)                      :: t_N
  real(kind=r_tran)                      :: t_S
  integer(kind=i_def), dimension(nfaces) :: opposite
  logical(kind=l_def), dimension(nfaces) :: missing_neighbour

  ! Zeroing of internal active variables
  dtdx = 0.0_r_tran
  dtdy = 0.0_r_tran
  t_N = 0.0_r_tran
  t_E = 0.0_r_tran
  t_S = 0.0_r_tran
  t_W = 0.0_r_tran

  ! For each face of cell, find the index in the neighbouring cell that
  ! corresponds to it.
  ! i.e for no orientation changes opposite = ( 3, 4, 1, 2 )
  ! We use the W2 map to determine these
  ! If there is no neighbour then we ensure the opposite points to
  ! the value on this edge
  opposite = -1
  missing_neighbour = .false.
  do df = 1, nfaces, 1
    df1 = map_w2(df)
    if (smap_w2_size(df) > 1) then
      ! There is a neighbour in direction df so find the
      ! neighboring edge corresponding to edge df
      do df2 = 1, nfaces, 1
        if (smap_w2(df2,2,df) == df1)  opposite(df) = df2
      end do
    else
      ! There is no neighbour in direction df so point to itself
      opposite(df) = df
      missing_neighbour(df) = .true.
    end if
  end do

  do k = nlayers - 1, 0, -1
    u =  0.5_r_tran*( wind(map_w2(1) + k) + wind(map_w2(3) + k) )
    v = -0.5_r_tran*( wind(map_w2(2) + k) + wind(map_w2(4) + k) )
    ik = 1 + k + (cell-1)*nlayers
    dtdx = dtdx + advective_increment(map_w3(1) + k) * u * real(m3_inv(ik,1,1), r_tran)
    dtdy = dtdy + advective_increment(map_w3(1) + k) * v * real(m3_inv(ik,1,1), r_tran)
    t_N = t_N + dtdy
    t_S = t_S - dtdy
    dtdy = 0.0_r_tran

    face = 4
    if ( v <= 0.0_r_tran .and. .not. missing_neighbour(face) ) then
      ! Use t_N from neighbouring column (if it exists)
      i_N = smap_md(1,2,face) + (opposite(face) - 1) * nlayers + k
    else
      ! Use t_N from this column
      i_N = map_md(1) + 3 * nlayers + k
    end if
    tracer(i_N) = tracer(i_N) + t_N
    t_N = 0.0_r_tran

    face = 2
    if ( v > 0.0_r_tran .and. .not. missing_neighbour(face) ) then
      ! Use t_S from neighbouring column (if it exists)
      i_S = smap_md(1,2,face) + (opposite(face) - 1) * nlayers + k
    else
      ! Use t_S from this column
      i_S = map_md(1) + nlayers + k
    end if
    tracer(i_S) = tracer(i_S) + t_S
    t_S = 0.0_r_tran

    t_E = t_E + dtdx
    t_W = t_W - dtdx
    dtdx = 0.0_r_tran

    face = 3
    if ( u <= 0.0_r_tran .and. .not. missing_neighbour(face) ) then
      ! Use t_E from neighbouring column (if it exists)
      i_E = smap_md(1,2,face) + (opposite(face) - 1) * nlayers + k
    else
      ! Use t_E from this column
      i_E = map_md(1) + 2 * nlayers + k
    end if
    tracer(i_E) = tracer(i_E) + t_E
    t_E = 0.0_r_tran

    face = 1
    if ( u > 0.0_r_tran .and. .not. missing_neighbour(face) ) then
      ! Use t_W from neighbouring column (if it exists)
      i_W = smap_md(1,2,face) + (opposite(face) - 1) * nlayers + k
    else
      ! Use t_W from this column
      i_W = map_md(1) + k
    end if
    tracer(i_W) = tracer(i_W) + t_W
    t_W = 0.0_r_tran
  end do

end subroutine adj_w3h_advective_update_code

end module adj_w3h_advective_update_kernel_mod
