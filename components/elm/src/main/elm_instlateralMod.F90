module elm_instlateralMod

  !-----------------------------------------------------------------------
  ! initialize elm data types
  !
#include <petsc/finclude/petsc.h>
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use abortutils   , only : endrun
  use decompMod    , only : bounds_type
  use GridCellType , only : grc_pp
  use elm_varpar   , only : nlevgrnd
  use petscsys
  use petscvec

  !-----------------------------------------
  ! Definition of component types
  !-----------------------------------------
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use ColumnType        , only : column_physical_properties
  use ColumnDataType    , only : column_water_state
  use GridcellType      , only : gridcell_physical_properties_type
  use elm_instMod       , only : soilstate_vars
  use elm_instMod       , only : col_pp
  use spmdMod           , only : masterproc, iam, npes, mpicom, comp_id

  implicit none
  save

  type(gridcell_physical_properties_type) :: ghost_grc_pp
  type(column_physical_properties)        :: ghost_col_pp
  type(column_water_state)                :: ghost_col_ws
  type(soilstate_type)                    :: ghost_soilstate_vars
  type(soilhydrology_type)                :: ghost_soilhydrology_vars

  public :: elm_instlateral_biophysics

contains

  
  !-----------------------------------------------------------------------
  subroutine elm_instlateral_biophysics(bounds_proc, domain_l)

    !
    ! DESCRIPTION
    ! initialize data structures
    use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
    use domainLateralMod     , only : domainlateral_type
    use UnstructuredGridType , only : ugdm_type, ugrid_type

    implicit none
    !
    type(bounds_type)                    :: bounds_proc
    type(domainlateral_type), intent(in) :: domain_l
    !
    type(bounds_type)                    :: bounds_ghost
    type(ugdm_type)  , pointer           :: ugdm
    type(ugrid_type) , pointer           :: ugrid
    PetscReal        , pointer           :: real_ptr(:)
    integer                              :: g,c,begc,endc

    ugdm => domain_l%dm_1dof
    ugrid => domain_l%ugrid

    call get_bounds_ghost(bounds_proc, bounds_ghost)

    call ghost_grc_pp%Init(bounds_proc%begg_ghost, bounds_proc%endg_ghost)

    call init_ghost_col_pp(bounds_proc, domain_l)

    allocate(ghost_col_ws%h2osoi_vol(bounds_proc%begc_ghost:bounds_proc%endc_ghost,1:nlevgrnd))
    ghost_col_ws%h2osoi_vol(:,:) = nan

    call ghost_soilhydrology_vars%InitAllocate(bounds_ghost)

    call init_ghost_soilstate_vars(bounds_proc, domain_l)

  end subroutine elm_instlateral_biophysics

  !------------------------------------------------------------------------
  subroutine get_bounds_ghost(bounds_proc, bounds_ghost)
    
    implicit none
    
    type(bounds_type), intent(in)  :: bounds_proc
    type(bounds_type), intent(out) :: bounds_ghost

    bounds_ghost%begg = bounds_proc%begg_ghost; bounds_ghost%endg = bounds_proc%endg_ghost
    bounds_ghost%begt = bounds_proc%begt_ghost; bounds_ghost%endt = bounds_proc%endt_ghost
    bounds_ghost%begl = bounds_proc%begl_ghost; bounds_ghost%endl = bounds_proc%endl_ghost
    bounds_ghost%begc = bounds_proc%begc_ghost; bounds_ghost%endc = bounds_proc%endc_ghost
    bounds_ghost%begp = bounds_proc%begp_ghost; bounds_ghost%endp = bounds_proc%endp_ghost

    bounds_ghost%begc_all = bounds_proc%begc_ghost; bounds_ghost%endc_all = bounds_proc%endc_ghost

  end subroutine get_bounds_ghost

  !------------------------------------------------------------------------
  subroutine init_ghost_col_pp(bounds_proc, domain_l)

    !
    ! DESCRIPTION
    ! Allocate memory for ghost cells and initialize by exchanging
    ! values from neigbhoring MPI ranks
    !
    use domainLateralMod     , only : domainlateral_type
    use UnstructuredGridType , only : ugdm_type, ugrid_type
    use UnstructuredGridType , only : ScatterDataG2L
    use spmdMod              , only : masterproc, iam, npes, mpicom, comp_id
    !
    implicit none
    !
    type(bounds_type)                    :: bounds_proc
    type(domainlateral_type), intent(in) :: domain_l
    !
    type(bounds_type)                    :: bounds_ghost
    type(ugdm_type)  , pointer           :: ugdm
    type(ugrid_type) , pointer           :: ugrid
    PetscInt                             :: g, c, j, idx
    PetscInt                             :: ndata_send, ndata_recv, nblocks
    PetscReal        , pointer           :: data_send(:), data_recv(:)

    ugdm => domain_l%dm_1dof
    ugrid => domain_l%ugrid

    call get_bounds_ghost(bounds_proc, bounds_ghost)
    call ghost_col_pp%Init(bounds_proc%begc_ghost, bounds_proc%endc_ghost)

    nblocks = nlevgrnd+1
    ndata_send = nblocks*domain_l%ugrid%ngrid_local
    ndata_recv = nblocks*domain_l%ugrid%ngrid_ghosted

    allocate(data_send(ndata_send))
    allocate(data_recv(ndata_recv))

    ! Note: The code below is exctracting data from the first "ngcells=endg-begg+1"
    !       soil columns. The assumpition is that each grid cell has a naturally-vegetated
    !       column
    idx = 0;
    do g = bounds_proc%begg, bounds_proc%endg
       c = bounds_proc%begc + g - bounds_proc%begg

       do j = 0, nlevgrnd
          idx = idx + 1; data_send(idx) = col_pp%zi(c,j)
       end do
    end do

    call ScatterDataG2L(domain_l%ugrid, nblocks, ndata_send, data_send, ndata_recv, data_recv)

    idx = (bounds_proc%endg - bounds_proc%begg + 1)*nblocks
    do g = bounds_ghost%begg, bounds_ghost%endg
       c = g
       do j = 0, nlevgrnd
          idx = idx + 1; ghost_col_pp%zi(c,j) = data_recv(idx)
       end do
    end do

    deallocate(data_send)
    deallocate(data_recv)

  end subroutine init_ghost_col_pp

  !------------------------------------------------------------------------
  subroutine init_ghost_soilstate_vars(bounds_proc, domain_l)

    !
    ! DESCRIPTION
    ! Allocate memory for ghost cells and initialize by exchanging
    ! values from neigbhoring MPI ranks
    !
    use domainLateralMod     , only : domainlateral_type
    use UnstructuredGridType , only : ugdm_type, ugrid_type
    use UnstructuredGridType , only : ScatterDataG2L
    use spmdMod              , only : masterproc, iam, npes, mpicom, comp_id
    !
    implicit none
    !
    type(bounds_type)                    :: bounds_proc
    type(domainlateral_type), intent(in) :: domain_l
    !
    type(bounds_type)          :: bounds_ghost
    type(ugdm_type)  , pointer :: ugdm
    type(ugrid_type) , pointer :: ugrid
    PetscReal        , pointer :: data_send(:), data_recv(:)
    integer                    :: g, c, j, idx
    PetscInt                   :: ndata_send, ndata_recv, nblocks

    ugdm => domain_l%dm_1dof
    ugrid => domain_l%ugrid

    call get_bounds_ghost(bounds_proc, bounds_ghost)

    nblocks = 3*nlevgrnd
    ndata_send = nblocks*domain_l%ugrid%ngrid_local
    ndata_recv = nblocks*domain_l%ugrid%ngrid_ghosted

    call ghost_soilstate_vars%InitAllocate(bounds_ghost)
    allocate(data_send(ndata_send))
    allocate(data_recv(ndata_recv))

    ! Note: The code below is exctracting data from the first "ngcells=endg-begg+1"
    !       soil columns. The assumpition is that each grid cell has a naturally-vegetated
    !       column
    idx = 0;
    do g = bounds_proc%begg, bounds_proc%endg
       c = bounds_proc%begc + g - bounds_proc%begg

       do j = 1, nlevgrnd
          idx = idx + 1; data_send(idx) = soilstate_vars%bsw_col(c,j)
          idx = idx + 1; data_send(idx) = soilstate_vars%watsat_col(c,j)
          idx = idx + 1; data_send(idx) = soilstate_vars%hksat_col(c,j)
       end do
    end do

    call ScatterDataG2L(domain_l%ugrid, nblocks, ndata_send, data_send, ndata_recv, data_recv)

    idx = (bounds_proc%endg - bounds_proc%begg + 1)*nblocks
    do g = bounds_ghost%begg, bounds_ghost%endg
       c = g
       do j = 1, nlevgrnd
          idx = idx + 1; ghost_soilstate_vars%bsw_col(c,j)    = data_recv(idx)
          idx = idx + 1; ghost_soilstate_vars%watsat_col(c,j) = data_recv(idx)
          idx = idx + 1; ghost_soilstate_vars%hksat_col(c,j)  = data_recv(idx)
       end do
    end do

    deallocate(data_send)
    deallocate(data_recv)

  end subroutine init_ghost_soilstate_vars

end module elm_instlateralMod
