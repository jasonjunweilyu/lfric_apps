!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!>
!> @brief   Set up specified mesh(es) from global/local mesh input file(s).
!> @details This routine will create a mesh_object_type(s) from a
!>          specified mesh input file and extrusion.
!>
!>          The algorithm differs depending on whether the input files(s)
!>          are prepartitioned (local mesh files) or not (global mesh files).
!>
!>          The result will be:
!>            * A set of local mesh objects stored in the application local
!>              mesh collection object.
!>            * A set of mesh objects stored in the application mesh collection
!>              object.
!>
!>          Local mesh object names will use the same mesh name as given in the
!>          input file. Extruded meshes are allowed to have an alternative mesh
!>          name as several meshes could exist in memory based on the same local
!>          mesh object.
!>
module lfric2lfric_init_mesh_mod

  use add_mesh_map_mod,            only: assign_mesh_maps
  use constants_mod,               only: str_def, i_def, l_def, str_max_filename
  use check_global_mesh_mod,       only: check_global_mesh
  use check_local_mesh_mod,        only: check_local_mesh
  use create_mesh_mod,             only: create_mesh
  use extrusion_mod,               only: extrusion_type
  use global_mesh_mod,             only: global_mesh_type
  use load_global_mesh_mod,        only: load_global_mesh
  use load_local_mesh_mod,         only: load_local_mesh
  use load_local_mesh_maps_mod,    only: load_local_mesh_maps
  use log_mod,                     only: log_event,         &
                                         log_scratch_space, &
                                         log_level_info,    &
                                         log_level_error,   &
                                         log_level_debug
  use namelist_collection_mod,     only: namelist_collection_type
  use namelist_mod,                only: namelist_type
  use partition_mod,               only: partitioner_interface
  use sci_query_mod,               only: check_uniform_partitions
  use runtime_partition_lfric_mod, only: get_partition_parameters
  use runtime_partition_mod,       only: mesh_cubedsphere,       &
                                         mesh_planar,            &
                                         create_local_mesh_maps, &
                                         create_local_mesh
  use global_mesh_collection_mod,  only: global_mesh_collection

  ! Configuration modules
  use finite_element_config_mod,   only: cellshape_quadrilateral
  use base_mesh_config_mod,        only: geometry_spherical, &
                                         topology_fully_periodic

  ! Lfric2lfric modules
  use lfric2lfric_config_mod,      only: regrid_method_map

  implicit none

  private
  public :: init_mesh

contains
!=======================================


!===============================================================================
!> @brief  Generates mesh(es) from mesh input file(s) on a given extrusion.
!>
!> @param[in] configuration          Application configuration object.
!>                                   This configuration object should contain the
!>                                   following defined namelist objects:
!>                                      * base_mesh
!>                                      * partititioning
!> @param[in] local_rank             The MPI rank of this process.
!> @param[in] total_ranks            Total number of MPI ranks in this job.
!> @param[in] destination_mesh_name  Destination mesh name to load from the mesh
!>                                   input file(s).
!> @param[in] source_mesh_name       Source mesh name to load from the mesh
!>                                   input file(s).
!> @param[in] extrusion              Extrusion object to be applied to meshes.
!> @param[in] stencil_depth          Required stencil depth for the application.
!> @param[in] regrid_method          Apply check for even partitions with the
!>                                   configured partition strategy if the
!>                                   regridding method is 'map'.
!>                                   (unpartitioned mesh input only)
!===============================================================================
subroutine init_mesh( configuration,           &
                      local_rank, total_ranks, &
                      destination_mesh_name,   &
                      source_mesh_name,        &
                      extrusion,               &
                      stencil_depth,           &
                      regrid_method )

  implicit none

  ! Arguments
  type(namelist_collection_type) :: configuration

  integer(kind=i_def),   intent(in) :: local_rank
  integer(kind=i_def),   intent(in) :: total_ranks
  character(len=*),      intent(in) :: destination_mesh_name
  character(len=*),      intent(in) :: source_mesh_name
  class(extrusion_type), intent(in) :: extrusion
  integer(kind=i_def),   intent(in) :: stencil_depth
  integer(kind=i_def),   intent(in) :: regrid_method

  ! Parameters
  character(len=9), parameter :: routine_name = 'init_mesh'

  ! Counters
  integer(kind=i_def) :: i

  ! Namelist variables
  type(namelist_type), pointer :: base_mesh_nml      => null()
  type(namelist_type), pointer :: partitioning_nml   => null()
  type(namelist_type), pointer :: lfric2lfric_nml    => null()

  logical(kind=l_def)     :: prepartitioned
  logical(l_def) :: generate_inner_haloes

  integer(kind=i_def) :: geometry
  integer(kind=i_def) :: topology
  integer(kind=i_def) :: mesh_selection

  character(len=str_max_filename)  :: destination_meshfile_prefix
  character(len=str_max_filename)  :: source_meshfile_prefix

  ! Local variables
  character(len=str_max_filename)     :: source_mesh_file
  character(len=str_max_filename)     :: destination_mesh_file
  character(len=str_def)              :: mesh_names(2)

  type(global_mesh_type),           pointer :: global_mesh_ptr => null()
  procedure(partitioner_interface), pointer :: partitioner_ptr => null()

  integer(kind=i_def) :: xproc  ! Processor ranks in mesh panel x-direction
  integer(kind=i_def) :: yproc  ! Processor ranks in mesh panel y-direction
  logical(kind=l_def) :: partitions_good

  !============================================================================
  ! Extract and check configuration variables
  !============================================================================
  base_mesh_nml      => configuration%get_namelist('base_mesh')
  partitioning_nml   => configuration%get_namelist('partitioning')
  lfric2lfric_nml    => configuration%get_namelist('lfric2lfric')

  call base_mesh_nml%get_value( 'prepartitioned', prepartitioned )
  call partitioning_nml%get_value( 'generate_inner_haloes', &
                                    generate_inner_haloes )
  call lfric2lfric_nml%get_value( 'destination_meshfile_prefix', &
                                   destination_meshfile_prefix )
  call lfric2lfric_nml%get_value( 'source_meshfile_prefix', &
                                   source_meshfile_prefix )

  if ( regrid_method == regrid_method_map .and. &
     trim(source_meshfile_prefix) /= trim(destination_meshfile_prefix) ) then

    write( log_scratch_space, '(A)' )                                &
         'When using LFRic intermesh maps, source and destination '//&
         'meshes should be extracted from the same file.'
    call log_event(log_scratch_space, log_level_error)
  end if

  mesh_names(1) = destination_mesh_name
  mesh_names(2) = source_mesh_name

  !===========================================================================
  ! Create local mesh objects:
  !  Two code pathes presented, either:
  !  1. The input files have been pre-partitioned.
  !     Meshes and are simply read from file and local mesh objects
  !     are populated.
  !  2. The input files have not been partitioned.
  !     Global meshes are loaded from file and partitioning is applied
  !     at runtime.  NOTE: This option is provided as legacy, and support
  !     is on a best endeavours basis.
  !===========================================================================
  if (prepartitioned) then

    !==========================================================================
    !  Read in local meshes / partition information / mesh maps
    !  direct from file.
    !==========================================================================
    !
    ! For this local rank, a mesh input file with a common base name
    ! of the following form should exist.
    !
    !   <input_basename>_<local_rank>_<total_ranks>.nc
    !
    ! Where 1 rank is assigned to each mesh partition.
    write(destination_mesh_file,'(A,2(I0,A))') &
        trim(destination_meshfile_prefix) // '_', local_rank, '-', &
                                           total_ranks, '.nc'

    write(source_mesh_file,'(A,2(I0,A))') &
        trim(source_meshfile_prefix) // '_', local_rank, '-',  &
                                           total_ranks, '.nc'

    ! Read in all local mesh data for this rank and
    ! initialise local mesh objects from them.
    !===========================================================
    ! Each partitioned mesh file will contain meshes of the
    ! same name as all other partitions.
    call log_event( 'Using pre-partitioned mesh file:', log_level_info )
    call log_event( '   '//trim(destination_mesh_file), log_level_info )
    call log_event( "Loading local mesh(es)", log_level_info )

    if (destination_mesh_file == source_mesh_file) then
      call load_local_mesh( destination_mesh_file, mesh_names )
    else
      call load_local_mesh( destination_mesh_file, destination_mesh_name )

      call log_event( 'Using pre-partitioned mesh file:', log_level_info )
      call log_event( '   '//trim(source_mesh_file), log_level_info )
      call log_event( "Loading local mesh(es)", log_level_info )

      call load_local_mesh( source_mesh_file, source_mesh_name )
    endif

    ! Apply configuration related checks to ensure that these
    ! meshes are suitable for the supplied application
    ! configuration.
    !===========================================================
    call check_local_mesh( configuration, &
                           stencil_depth, &
                           mesh_names )

    ! Load and assign mesh maps.
    !===========================================================
    ! Mesh map identifiers are determined by the source/target
    ! mesh IDs they relate to. As a result inter-grid mesh maps
    ! need to be loaded after the relevant local meshes have
    ! been loaded.
    if (regrid_method == regrid_method_map) then
      call load_local_mesh_maps( destination_mesh_file, mesh_names )
    end if
  else

    !==========================================================================
    ! Perform runtime partitioning of global meshes.
    !==========================================================================
    call base_mesh_nml%get_value( 'geometry', geometry )
    call base_mesh_nml%get_value( 'topology', topology )

    if ( geometry == geometry_spherical .and. &
         topology == topology_fully_periodic ) then
      mesh_selection = mesh_cubedsphere
      call log_event( "Setting up cubed-sphere partition mesh(es)", &
                      log_level_debug )
    else
      mesh_selection = mesh_planar
      call log_event( "Setting up planar partition mesh(es)", &
                      log_level_debug )
    end if

    call log_event( "Setting up partition mesh(es)", log_level_info )
    write(source_mesh_file,'(A)') trim(source_meshfile_prefix) // '.nc'
    write(destination_mesh_file,'(A)') trim(destination_meshfile_prefix) &
                                                               // '.nc'

    ! Set constants that will control partitioning.
    !===========================================================
    call get_partition_parameters( configuration, mesh_selection, &
                                   total_ranks, xproc, yproc,     &
                                   partitioner_ptr )

    ! Read in all global meshes from input file
    !===========================================================
    if (destination_mesh_file == source_mesh_file) then
      call load_global_mesh( destination_mesh_file, mesh_names )
    else
      call load_global_mesh( destination_mesh_file, destination_mesh_name )
      call load_global_mesh( source_mesh_file, destination_mesh_name )
    endif

    if (regrid_method == regrid_method_map) then

      ! Apply configuration related checks to ensure that these
      ! meshes are suitable for the supplied application
      ! configuration.
      !===========================================================
      call check_global_mesh( configuration, mesh_names )

      ! Optional, Check for even partitions
      !===========================================================
      ! Check global meshes which may require even partition sizes,
      ! these will generally be where there are meshes of differing
      ! resolutions and the node locations of one mesh overlay the
      ! the other(s). Checking the lowest resolution mesh for
      ! even partitions should ensure common geographical partitions.
      do i=1, size(mesh_names)

        write(log_scratch_space, '(A)') &
            'Checking partitioning on mesh '//trim(mesh_names(i))
        call log_event( log_scratch_space, log_level_info )

        global_mesh_ptr => &
                global_mesh_collection%get_global_mesh( mesh_names(i) )
        partitions_good = check_uniform_partitions( global_mesh_ptr, &
                                                    xproc, yproc )
        if ( .not. partitions_good ) then
          write(log_scratch_space, '(A)')                                &
              'Global mesh ('// trim(global_mesh_ptr%get_mesh_name()) // &
              ') ' // 'does not produce even partitions with the ' //    &
              'requested partition strategy'
          call log_event( log_scratch_space, log_level_error )
        end if
      end do
    end if

    ! Partition the global meshes
    !===========================================================
    call create_local_mesh( mesh_names,              &
                            local_rank, total_ranks, &
                            xproc, yproc,            &
                            stencil_depth,           &
                            generate_inner_haloes,   &
                            partitioner_ptr )

    ! Read in the global intergrid mesh mappings,
    ! then create the associated local mesh maps
    !===========================================================
    if (regrid_method == regrid_method_map) then
      call create_local_mesh_maps( destination_mesh_file )
    end if

  end if  ! prepartitioned

  !============================================================================
  ! Extrude the specified meshes from local mesh objects into
  ! mesh objects on the given extrusion.
  ! Alternative names are needed in case the source and destination
  ! mesh files use the same mesh name.
  !============================================================================
  call create_mesh( mesh_names, extrusion )

  !============================================================================
  ! Generate intergrid LiD-LiD maps and assign them to mesh objects.
  !============================================================================
  if (regrid_method == regrid_method_map) then
    call assign_mesh_maps(mesh_names)
  end if

end subroutine init_mesh

end module lfric2lfric_init_mesh_mod
