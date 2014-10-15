module VariableSet_mod

  ! モジュール引用; Use statements
  !

  use dc_types, only: DP

  use dc_message, only: MessageNotify

  use GeometricField_mod

  use GridSet_mod, only: plMesh

  ! 宣言文; Declareration statements
  !

  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: VariableSet_Init, VariableSet_Final

  type(volScalarField), public :: v_h, v_U1, v_U2, v_BtmTopl
  real(DP), public, dimension(:,:), allocatable :: wc_h, wc_hU1, wc_hU2
  real(DP), public, dimension(:,:), allocatable :: wc_U1, wc_U2, wc_BtmTopl, wc_dhdt, wc_etc
  real(DP), public, dimension(:,:), allocatable :: wc_WindStress1, wc_WindStress2

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'VariableSet_mod' !< Module Name

contains

subroutine VariableSet_Init()

  !
  use DGHelper_mod, only: DGHelper_MallocElemNode


  !

  call DGHelper_MallocElemNode(wc_h)
  call DGHelper_MallocElemNode(wc_BtmTopl)
  call DGHelper_MallocElemNode(wc_hU1)
  call DGHelper_MallocElemNode(wc_U1)
  call DGHelper_MallocElemNode(wc_hU2)
  call DGHelper_MallocElemNode(wc_U2)
  call DGHelper_MallocElemNode(wc_dhdt)
  call DGHelper_MallocElemNode(wc_etc)
  
  call DGHelper_MallocElemNode(wc_WindStress1)
  call DGHelper_MallocElemNode(wc_WindStress2)

  contains

end subroutine VariableSet_Init

subroutine VariableSet_Final()

  if(allocated(wc_h)) then
     deallocate(wc_h, wc_BtmTopl, wc_hU1, wc_hU2)
     deallocate(wc_U1, wc_U2, wc_dhdt)
     deallocate(wc_WindStress1, wc_WindStress2)
  end if

end subroutine VariableSet_Final


end module VariableSet_mod
