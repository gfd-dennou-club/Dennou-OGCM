!-------------------------------------------------------------
! Copyright (c) 2013-2014 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module BudgetAnalysis_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING, TOKEN

  use dc_message, only: &
       & MessageNotify
  
  use gtool_history

  use Constants_mod, only: &
       & RPlanet, Grav, RefDens, &
       & hViscCoef, vViscCoef, hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use HydroBouEqSolverRHS_mod
  use HydroBouEqSolverVImplProc_mod

  use DiagVarFileSet_mod
  use DiagnoseUtil_mod
  use DiagVarEval_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: BudgetAnalysis_Init, BudgetAnalysis_Final
  public :: BudgetAnalysis_perform

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'BudgetAnalysis_mod' !< Module Name

  character(*), parameter, public :: BUDGETANAKEY_GLOBALMEANENERGY = 'GlobalMeanEnergy'
  character(*), parameter, public :: BUDGETANAKEY_ENERGYBUDGET = 'EnergyBudget'
  character(*), parameter, public :: BUDGETANAKEY_TEAVG = 'TEAvg'
  character(*), parameter, public :: BUDGETANAKEY_KEAVG = 'KEAvg'
  character(*), parameter, public :: BUDGETANAKEY_PEAVG = 'PEAvg'
  character(*), parameter, public :: BUDGETANAKEY_PE2KEAVG = 'PE2KEAvg'
  character(*), parameter, public :: BUDGETANAKEY_KE2PEAVG = 'KE2PEAvg'
  character(*), parameter, public :: BUDGETANAKEY_KEINPUTSURFAVG = 'KEInputSurfAvg'
  character(*), parameter, public :: BUDGETANAKEY_HDIFFDISPAVG = 'HDiffDispAvg'
  character(*), parameter, public :: BUDGETANAKEY_VDIFFDISPAVG = 'VDiffDispAvg'
  character(*), parameter, public :: BUDGETANAKEY_ADVECTWORKAVG = 'AdvectWorkAvg'
  character(*), parameter, public :: BUDGETANAKEY_KEGENNETAVG = 'KEGenNet'


  character(*), parameter, public :: BUDGETANAKEY_ANGMOMBUDGET = 'AngMomBudget'
  character(*), parameter, public :: BUDGETANAKEY_ANGMOMAVG = 'AngMomAvg'

  type(gt_history) :: hst_globalMeanEnergy
  logical :: globalMeanEnergyFlag

  type(gt_history) :: hst_energyBudget
  logical :: energyBudgAnaFlag

  type(gt_history) :: hst_angMomBudget
  logical :: angMomBudgAnaFlag

contains

  !>
  !!
  !!
  subroutine BudgetAnalysis_Init(diagVar_gthsInfo, budgetAnaName)

    ! 宣言文; Declare statements
    !
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo
    character(*), intent(in) :: budgetAnaName(:)
    
    !
    !
    integer :: n

    ! 実行文; Executable statements
    !

    globalMeanEnergyFlag = .false.
    energyBudgAnaFlag = .false.
    angMomBudgAnaFlag = .false.

    do n=1, size(budgetAnaName)
       select case(budgetAnaName(n))
       case (BUDGETANAKEY_GLOBALMEANENERGY)
          globalMeanEnergyFlag = .true.
       case (BUDGETANAKEY_ENERGYBUDGET) 
          energyBudgAnaFlag = .true.
       case (BUDGETANAKEY_ANGMOMBUDGET) 
          angMomBudgAnaFlag = .true.
       case ('')
       case Default
          call MessageNotify('E', module_name, &
               & "The specified type of budget analysis '%c' is invalid.", c1=trim(budgetAnaName(n)) )
       end select
    end do

    call prepair_Output(diagVar_gthsInfo)

    call HydroBouEqSolverRHS_Init()
    call HydroBouEqSolverVImplProc_Init()

  end subroutine BudgetAnalysis_Init

  !>
  !!
  !!
  subroutine BudgetAnalysis_Final()

    ! 実行文; Executable statements
    !
    if( globalMeanEnergyFlag ) then
       call HistoryClose(hst_globalMeanEnergy)
    end if

    if( energyBudgAnaFlag ) then
       call HistoryClose(hst_energyBudget)
    end if

    if( angMomBudgAnaFlag ) then
       call HistoryClose(hst_angMomBudget)
    end if

    call HydroBouEqSolverRHS_Final()
    call HydroBouEqSolverVImplProc_Final()

  end subroutine BudgetAnalysis_Final

  subroutine BudgetAnalysis_perform( )

    use VariableSet_mod, only: &
         & xyz_UN, xyz_VN, xyz_PTempEddN, xyz_PTempEddB, xy_SurfHeightN, &
         & xy_totDepthBasic

    ! 宣言文; Declare statements
    !
 
    ! 実行文; Executable statements
    !
    
    if( globalMeanEnergyFlag) call analyze_globalMeanEnergy()
    if( energyBudgAnaFlag ) call analyze_energyBudget()
    if( angMomBudgAnaFlag ) call analyze_angMomBudget()

  contains
    subroutine analyze_globalMeanEnergy()

      use DiagVarSet_mod, only: xyz_DensEdd

      real(DP) :: KEnAvg, PEnAvg
      real(DP) :: xy_totDepth(0:iMax-1,jMax)

      call MessageNotify("M", module_name, "Calculate global mean of each energy ..")

      xy_totDepth = xy_totDepthBasic + xy_SurfHeightN

      KEnAvg = eval_kineticEnergyAvg(xyz_UN, xyz_VN)
      PEnAvg = eval_potentialEnergyAvg(xyz_DensEdd, xy_totDepth) 
      call HistoryPut( BUDGETANAKEY_KEAVG, KEnAvg, hst_globalMeanEnergy )
      call HistoryPut( BUDGETANAKEY_PEAVG, PEnAvg, hst_globalMeanEnergy )
      call HistoryPut( BUDGETANAKEY_TEAVG, KEnAvg + PEnAvg, hst_globalMeanEnergy )
    end subroutine analyze_globalMeanEnergy

    subroutine analyze_energyBudget()

use DiagVarSet_mod, only: xyz_DensEdd

use HydroBouEqSolverRHS_mod

      real(DP) :: KEnAvg, PEnAvg
      
      real(DP), dimension(2) :: PEConvert, HDiffDisp, VDiffDisp, KEInput, Advect
      real(DP) :: wz_Zero(lMax,0:kMax)
      real(DP), dimension(0:iMax-1,jMax,0:kMax,2) :: xyz_dudt, xyz_dvdt, xyz_dPTempdt
      real(DP) :: xy_totDepth(0:iMax-1,jMax)


      call MessageNotify("M", module_name, "Analyze energy budget..")
      xy_totDepth = xy_totDepthBasic

      !
      KEnAvg = eval_kineticEnergyAvg(xyz_UN, xyz_VN)
      PEnAvg = eval_potentialEnergyAvg(xyz_DensEdd, xy_totDepth) 
      call HistoryPut( BUDGETANAKEY_KEAVG, KEnAvg, hst_energyBudget )
      call HistoryPut( BUDGETANAKEY_PEAVG, PEnAvg, hst_energyBudget )
      call HistoryPut( BUDGETANAKEY_TEAVG, KEnAvg + PEnAvg, hst_energyBudget )


      KEInput=0d0; PEConvert=0d0; HDiffDisp=0d0; VDiffDisp=0d0; Advect = 0d0;
      call calc_dudtdvdt(xyz_dudt(:,:,:,1), xyz_dvdt(:,:,:,1), xyz_dPTempdt(:,:,:,1), &
           & KEInput(1), PEConvert(1), HDiffDisp(1), VDiffDisp(1), Advect(1))

!!
      call HistoryPut( BUDGETANAKEY_PE2KEAVG, PEConvert(1), hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_KE2PEAVG, -PEConvert(1), hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_KEINPUTSURFAVG, KEInput(1), hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_HDIFFDISPAVG, HDiffDisp(1), hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_VDIFFDISPAVG, VDiffDisp(1), hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_ADVECTWORKAVG, Advect(1), hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_KEGENNETAVG, & 
           & sum(PEConvert + KEInput + VDiffDisp + HDiffDisp + Advect), &
           & hst_energyBudget)

write(*,*) "HDiffDisp: VDiffDisp: PE2KE:  KEInput: Advect:"
write(*,'(5(E15.5e4,a))'), HDiffDisp(1), ":", VDiffDisp(1), ":", PEConvert(1), ":", KEInput(1), ":", Advect(1)
write(*,*) ":  KEGenNet", PEConvert(1)+KEInput(1)+VDiffDisp(1)+HDiffDisp(1)
write(*,*) ":  KEGenNet", PEConvert(1)+KEInput(1)+VDiffDisp(1)+HDiffDisp(1)+Advect(1)

    end subroutine analyze_energyBudget

    subroutine analyze_angMomBudget()

      real(DP) :: AngMomAvg

      call MessageNotify("M", module_name, "Analyze angular momentum budget..")

      AngMomAvg = eval_angularMomAvg(xyz_UN)
      call HistoryPut( BUDGETANAKEY_ANGMOMAVG, AngMomAvg, hst_angMomBudget )

    end subroutine analyze_angMomBudget

    subroutine calc_dudtdvdt(xyz_dudt, xyz_dvdt, xyz_dPTempdt, KEInput, PEConvert, HDiffDisp, VDiffDisp, Advect )

!use wa_module
use HydroBouEqSolverRHS_mod
use VariableSet_mod!, only: z_PTempBasic, xy_WindStressU, xy_WindStressV
use TemporalIntegSet_mod, only: DelTime
use BoundCondSet_mod, only: &
       & KinBC_Surface, DynBC_Surface, DynBC_Bottom, &
       & DynBCTYPE_Slip, DynBCTYPE_NoSlip
!use HydroBouEqSolver_mod

      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: xyz_dudt, xyz_dvdt, xyz_dPTempdt
      real(DP), intent(inout) :: KEInput, PEConvert, HDiffDisp, VDiffDisp, Advect
      

      !
      real(DP) :: xyz_DuDSig(0:iMax-1,jMax,0:kMax), xyz_DvDSig(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_Div(0:iMax-1,jMax,0:kMax), xyz_Vor(0:iMax-1,jMax,0:kMax)
      real(DP) , dimension(0:iMax-1,jMax,0:kMax) :: xyz_InvisRHSU, xyz_InvisRHSV, &
           & xyz_HDiffU, xyz_HDiffV, xyz_VDiff, xyz_PressEddTmp, &
           & xyz_UTmp, xyz_VTmp, xyz_PressEdd, xyz_PressWork, xyz_AdvWork, &
           & xyz_PTempBasic, xyz_DensEdd, xyz_dKdt
      real(DP) :: xyz_CosLat(0:iMax-1,jMax,0:kMax)
      real(DP), dimension(lMax, 0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_TmpChi, wz_TmpPsi, wz_Tmp
      real(DP), dimension(lMax, 0:kMax) :: wz_VorTmp, wz_DivTmp, wz_PTempEddTmp, wt_Vor, wt_Div, wt_PTempEdd, wz_VorA, wz_DivA
      real(DP) :: xy_totDepth(0:iMax-1,jMax)
      
      !
      xyz_CosLat = cos(xyz_Lat)
      xy_totDepth = xy_totDepthBasic

      xyz_Div = eval_Div(xyz_UN, xyz_VN)
      xyz_Vor = eval_Vor(xyz_UN, xyz_VN)

      !
      xyz_HDiffU = hDiffCoef * xyz_CosLat*xyz_AlphaOptr_wz(wz_xyz(xyz_Div), -wz_xyz(xyz_Vor))
      xyz_HDiffV = hDiffCoef * xyz_CosLat*xyz_AlphaOptr_wz(wz_xyz(xyz_Vor),  wz_xyz(xyz_Div))

      !
      xyz_SaltN = 0d0
      xy_totDepth = xy_totDepthBasic
      xyz_PTempBasic = spread(spread(z_PTempBasic,1,jMax),1,iMax)

      ! 
      ! Step 1
      xyz_DensEdd = eval_DensEdd(xyz_PTempBasic + xyz_PTempEddN, xyz_SaltN, xy_totDepth)
      xyz_SigDot = diagnose_SigDot( xy_totDepth, xyz_UN*xyz_CosLat, xyz_VN*xyz_CosLat, xyz_Div )
      wz_Tmp = 0d0
      xyz_PressEddTmp = Diagnose_HydroPressEdd(xy_totDepth, xyz_DensEdd)
      xyz_PressEdd = 0.5d0*xyz_PressEddTmp

!!$      call calc_VorEqDivEqInvisRHS(wz_VorRHS, wz_DivRHS, &
!!$         & xyz_Vor, xyz_UN*xyz_CosLat, xyz_VN*xyz_CosLat, xy_w(wz_Tmp(:,0)), xyz_DensEdd, &
!!$         & xyz_PressEddTmp, Diagnose_GeoPot(xy_totDepth), xyz_SigDot)

!!$      call correct_DivEqRHSUnderRigidLid(wz_DivRHS, &
!!$           & xy_SurfPressN, wz_xyz(xyz_Div), xy_totDepth, vDiffCoef, DelTime)

      wz_TmpPsi = wz_InvLapla2D_wz( wz_VorRHS  )
      wz_TmpChi = wz_InvLapla2D_wz( wz_DivRHS  ) + wz_xyz(xyz_PressEddTmp/RefDens)
      xyz_InvisRHSU = 0.5d0*xyz_CosLat*xyz_AlphaOptr_wz(wz_TmpChi, -wz_TmpPsi)
      xyz_InvisRHSV = 0.5d0*xyz_CosLat*xyz_AlphaOptr_wz(wz_TmpPsi,  wz_TmpChi)

      call calc_TracerEqInvisRHS( wz_PTempRHS, &
           & xyz_PTempBasic+xyz_PTempEddN, xyz_UN*xyz_CosLat, xyz_VN*xyz_CosLat, xyz_Div, xyz_SigDot )

      wz_VorTmp = wz_xyz( 4d0/12d0*eval_Vor(xyz_UB, xyz_VB) + 8d0/12d0*xyz_Vor ) + 10d0/12d0*DelTime*wz_VorRHS
      wz_DivTmp = wz_xyz( 4d0/12d0*eval_Div(xyz_UB, xyz_VB) + 8d0/12d0*xyz_Div ) + 10d0/12d0*DelTime*wz_DivRHS
      wz_PTempEddTmp = wz_xyz( 4d0/12d0*xyz_PTempEddB + 8d0/12d0*xyz_PTempEddN ) + 10d0/12d0*DelTime*wz_PTempRHS
    

      wt_Vor = wt_wz(wz_VorTmp); wt_Div = wt_wz(wz_DivTmp); wt_PTempEdd = wt_wz(wz_PTempEddTmp)
!      call apply_boundaryConditions(wt_Vor, wt_Div, wt_PTempEdd)
      wz_VorTmp = wz_wt(wt_Vor); wz_DivTmp = wz_wt(wt_Div); wz_PTempEddTmp = wz_wt(wt_PTempEdd)

      wz_TmpPsi = wz_InvLapla2D_wz( wz_VorTmp  )
      wz_TmpChi = wz_InvLapla2D_wz( wz_DivTmp  )
      xyz_UTmp =xyz_CosLat*xyz_AlphaOptr_wz(wz_TmpChi, -wz_TmpPsi)
      xyz_VTmp =xyz_CosLat*xyz_AlphaOptr_wz(wz_TmpPsi,  wz_TmpChi)

      ! Step 2
      !
      xyz_DensEdd = eval_DensEdd(xyz_PTempBasic + xyz_wz(wz_PTempEddTmp), xyz_SaltN, xy_totDepth)
      xyz_SigDot = diagnose_SigDot( xy_totDepth, xyz_UTmp*xyz_CosLat, xyz_VTmp*xyz_CosLat, xyz_wz(wz_DivTmp) )
      xyz_PressEddTmp = Diagnose_HydroPressEdd(xy_totDepth, xyz_DensEdd)
      xyz_PressEdd = xyz_PressEddTmp! + 0.5d0*xyz_PressEddTmp
      wz_Tmp = 0d0
!!$      call calc_VorEqDivEqInvisRHS(wz_VorRHS, wz_DivRHS, &
!!$         & xyz_wz(wz_VorTmp), xyz_UTmp*xyz_CosLat, xyz_VTmp*xyz_CosLat, xy_w(wz_Tmp(:,0)), xyz_DensEdd, &
!!$         & xyz_PressEddTmp, Diagnose_GeoPot(xy_totDepth), xyz_SigDot)

      wz_TmpPsi = wz_InvLapla2D_wz( wz_VorRHS  )
      wz_TmpChi = wz_InvLapla2D_wz( wz_DivRHS  ) + wz_xyz(xyz_PressEddTmp/RefDens)
      xyz_InvisRHSU = xyz_CosLat*xyz_AlphaOptr_wz(wz_TmpChi, -wz_TmpPsi)
      xyz_InvisRHSV = xyz_CosLat*xyz_AlphaOptr_wz(wz_TmpPsi,  wz_TmpChi)

!!$      call calc_HydroBouEqHViscRHS(wz_VorRHS, wz_DivRHS, wz_PTempRHS, &
!!$           & wz_VorTmp, wz_DivTmp, wz_PTempEddTmp, &
!!$           & hDiffCoef, hHyperViscCoef, hDiffCoef, &
!!$           & isRHSReplace=.false. )

!!$      call calc_HydroBouEqVViscRHS(wz_VorRHS, wz_DivRHS, wz_PTempRHS, &
!!$           & wz_VorTmp, wz_DivTmp, wz_PTempEddTmp, &
!!$           & 0.5d0*vDiffCoef, vHyperViscCoef, 0.5d0*vDiffCoef, &
!!$           & isRHSReplace=.false. )

!!$      call correct_DivEqRHSUnderRigidLid(wz_DivRHS, &
!!$           & xy_SurfPressN, wz_DivTmp, xy_totDepth, vDiffCoef, DelTime)

      wz_VorA = wz_xyz(xyz_Vor) + wz_VorRHS*DelTime
      wz_DivA = wz_xyz(xyz_Div) + wz_DivRHS*DelTime
!!$      wz_VorA = wz_VorA + wz_VorRHS*0.5d0*DelTime
!!$      wz_DivA = wz_DivA + wz_DivRHS*0.5d0*DelTime

      wt_Vor = wt_wz(wz_VorA); wt_Div = wt_wz(wz_DivA); 
!!$      call Advance_VImplicitProc( wt_Vor, wt_Div, wt_PTempEdd, &
!!$           & xy_WindStressU, xy_WindStressV, xy_totDepth, &
!!$           & 0.5d0, DelTime, &
!!$           & DynBC_Surface, DynBC_Bottom )

      wz_TmpPsi = wz_InvLapla2D_wz( wz_wt(wt_Vor) )
      wz_TmpChi = wz_InvLapla2D_wz( wz_wt(wt_Div) )
      xyz_UA =xyz_CosLat*xyz_AlphaOptr_wz(wz_TmpChi, -wz_TmpPsi)
      xyz_VA =xyz_CosLat*xyz_AlphaOptr_wz(wz_TmpPsi,  wz_TmpChi)

      !

      xyz_UTmp = 0.5d0*(xyz_UN + xyz_UA); xyz_VTmp = 0.5d0*(xyz_VN + xyz_VA)

      HDiffDisp = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
           & xyz_UTmp*xyz_HDiffU + xyz_VTmp*xyz_HDiffV &
           & ))

      wz_Tmp = 0d0
      xyz_PressEdd = xyz_PressEdd + spread(xy_SurfPressN,3,kMax+1)
      xyz_PressWork = - ( &
           &   xyz_UTmp * xyz_CosLat*xyz_AlphaOptr_wz(wz_xyz(xyz_PressEdd/RefDens), wz_Tmp) &
           & + xyz_VTmp * xyz_CosLat*xyz_AlphaOptr_wz(wz_Tmp, wz_xyz(xyz_PressEdd/RefDens)) &
           & )

      PEConvert = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyz_PressWork ))


      xyz_AdvWork = xyz_UTmp*xyz_InvisRHSU + xyz_VTmp*xyz_InvisRHSV 
      Advect = AvrLonLat_xy(xy_IntSig_BtmToTop_xyz( xyz_AdvWork ))


!!$      xyz_VDiff = vDiffCoef*( &
!!$           &    xyz_UTmp * xyz_wt(wt_DSig_wt(wt_DSig_wt(wt_xyz(xyz_UA)))) &
!!$           &  + xyz_VTmp * xyz_wt(wt_DSig_wt(wt_DSig_wt(wt_xyz(xyz_VA)))) &
!!$           & )/spread(xy_totDepth**2,3,kMax+1)

      xyz_VDiff = - vDiffCoef*( &
           &    xyz_wt(wt_DSig_wt(wt_xyz(xyz_UTmp)))**2 &
           &  + xyz_wt(wt_DSig_wt(wt_xyz(xyz_VTmp)))**2 &
           & )/spread(xy_totDepth**2,3,kMax+1)

      VDiffDisp =  AvrLonLat_xy( xy_IntSig_BtmToTop_xyz(xyz_VDiff) )
!      xyz_VDiff(:,:,0) = 0.5d0*( (xyz_UA(:,:,0)**2+xyz_VA(:,:,0)**2) - (xyz_UN(:,:,0)**2+xyz_VN(:,:,0)**2))/DelTime


      xyz_dKdt = 0.5d0*((xyz_UA + xyz_UN)*(xyz_UA - xyz_UN) + (xyz_VA + xyz_VN)*(xyz_VA - xyz_VN))/DelTime
      KEInput = AvrLonLat_xy(xy_IntSig_BtmToTop_xyz( xyz_dKdt )) - PEConvert - Advect - VDiffDisp - HDiffDisp

write(*,*) "*KEInut surf:", KEInput

write(*,*) "* DKEDt=", AvrLonLat_xy(xy_IntSig_BtmToTop_xyz( xyz_dKdt ))

    end subroutine calc_dudtdvdt

  end subroutine BudgetAnalysis_perform

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine prepair_Output(diagVar_gthsInfo)
    ! 宣言文; Declaration statement
    !
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    if ( globalMeanEnergyFlag ) then
       call HistoryCreate( & 
            & file= trim(diagVar_gthsInfo%FilePrefix) // 'GlobalMeanEnergy.nc', title='global mean of each energy', &
            & source='Dennou-OGCM', &
            & institution='Dennou-OGCM project', &
            & dims=(/'t'/), dimsizes=(/ 0 /), &
            & longnames=(/'time'/),&
            & units=(/ diagVar_gthsInfo%intUnit /), &
            & origin=real(diagVar_gthsInfo%origin), interval=real(diagVar_gthsInfo%intValue), &
            & history=hst_globalMeanEnergy )  

write(*,*) diagVar_gthsInfo%origin

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_TEAVG, dims=(/'t'/), &
             & longname='global mean of total energy', units='J*m-3*kg-1', xtype='double',&
             & history=hst_globalMeanEnergy)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_KEAVG, dims=(/'t'/), &
             & longname='global mean of kinetic energy', units='J*m-3*kg-1', xtype='double',&
             & history=hst_globalMeanEnergy)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_PEAVG, dims=(/'t'/), &
             & longname='global mean of pertubated potential energy', units='J*m-3*kg-1', xtype='double',&
             & history=hst_globalMeanEnergy)
     end if

    if ( energyBudgAnaFlag ) then
       call HistoryCreate( & 
            & file= trim(diagVar_gthsInfo%FilePrefix) // 'EnergyBudget.nc', title='energy budget analysis', &
            & source='Dennou-OGCM', &
            & institution='Dennou-OGCM project', &
            & dims=(/'t'/), dimsizes=(/ 0 /), &
            & longnames=(/'time'/),&
            & units=(/ diagVar_gthsInfo%intUnit /), &
            & origin=real(diagVar_gthsInfo%origin), interval=real(diagVar_gthsInfo%intValue), &
            & history=hst_energyBudget )  

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_TEAVG, dims=(/'t'/), &
             & longname='global mean of total energy', units='J*m-3*kg-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_KEAVG, dims=(/'t'/), &
             & longname='global mean of kinetic energy', units='J*m-3*kg-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_PEAVG, dims=(/'t'/), &
             & longname='global mean of pertubated potential energy', units='J*m-3*kg-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_PE2KEAVG, dims=(/'t'/), &
             & longname='global mean of conversion of potential energy into kinetic energy', &
             & units='J*m-3*kg-1*s-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_KE2PEAVG, dims=(/'t'/), &
             & longname='global mean of conversion of potential energy into kinetic energy', &
             & units='J*m-3*kg-1*s-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_KEINPUTSURFAVG, dims=(/'t'/), &
             & longname='global mean of the input of kinetic energy on the sea surface', &
             & units='J*m-3*kg-1*s-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_HDIFFDISPAVG, dims=(/'t'/), &
             & longname='global mean of the dissipation due to horizontal diffusion', &
             & units='J*m-3*kg-1*s-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_VDIFFDISPAVG, dims=(/'t'/), &
             & longname='global mean of the dissipation due to vertical diffusion', &
             & units='J*m-3*kg-1*s-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_ADVECTWORKAVG, dims=(/'t'/), &
             & longname='global mean of the kinetic energy generation due to advection', &
             & units='J*m-3*kg-1*s-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_KEGENNETAVG, dims=(/'t'/), &
             & longname='global mean of net kinetic energy generation', &
             & units='J*m-3*kg-1*s-1', xtype='double',&
             & history=hst_energyBudget)
       
    end if

    if ( angMomBudgAnaFlag ) then

       call HistoryCreate( & 
            & file= trim(diagVar_gthsInfo%FilePrefix) // 'AngMomBudget.nc', title='energy budget analysis', &
            & source='Dennou-OGCM', &
            & institution='Dennou-OGCM project', &
            & dims=(/'t'/), dimsizes=(/ 0 /), &
            & longnames=(/'time'/),&
            & units=(/ diagVar_gthsInfo%intUnit /), &
            & origin=real(0), interval=real(diagVar_gthsInfo%intValue), &
            & history=hst_angMomBudget )  

       call HistoryAddVariable( &
            & varname=BUDGETANAKEY_ANGMOMAVG, dims=(/'t'/), & 
            & longname='global mean of relative angular momentum', units='kg*m2*s-1', xtype='double',&
             & history=hst_angMomBudget)
    end if

  end subroutine prepair_Output

end module BudgetAnalysis_mod

