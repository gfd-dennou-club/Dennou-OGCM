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
       & hDiffCoef, vDiffCoef

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

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

  character(*), parameter, public :: BUDGETANAKEY_ENERGYBUDGET = 'EnergyBudget'
  character(*), parameter, public :: BUDGETANAKEY_TEAVG = 'TEAvg'
  character(*), parameter, public :: BUDGETANAKEY_KEAVG = 'KEAvg'
  character(*), parameter, public :: BUDGETANAKEY_PEAVG = 'PEAvg'
  character(*), parameter, public :: BUDGETANAKEY_PE2KEAVG = 'PE2KEAvg'
  character(*), parameter, public :: BUDGETANAKEY_KE2PEAVG = 'KE2PEAvg'
  character(*), parameter, public :: BUDGETANAKEY_KEINPUTSURFAVG = 'KEInputSurfAvg'
  character(*), parameter, public :: BUDGETANAKEY_HDIFFDISPAVG = 'HDiffDispAvg'
  character(*), parameter, public :: BUDGETANAKEY_VDIFFDISPAVG = 'VDiffDispAvg'
  character(*), parameter, public :: BUDGETANAKEY_SURFPRESSWORKAvg = 'SurfPressWorkAvg'
  character(*), parameter, public :: BUDGETANAKEY_KEGENNETAVG = 'KEGenNet'


  character(*), parameter, public :: BUDGETANAKEY_ANGMOMBUDGET = 'AngMomBudget'
  character(*), parameter, public :: BUDGETANAKEY_ANGMOMAVG = 'AngMomAvg'

  logical :: energyBudgAnaFlag
  logical :: angMomBudgAnaFlag
  type(gt_history) :: hst_energyBudget
  type(gt_history) :: hst_angMomBudget

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

    energyBudgAnaFlag = .false.
    angMomBudgAnaFlag = .false.

    do n=1, size(budgetAnaName)
       select case(budgetAnaName(n))
       case (BUDGETANAKEY_ENERGYBUDGET) 
          energyBudgAnaFlag = .true.
       case (BUDGETANAKEY_ANGMOMBUDGET) 
          angMomBudgAnaFlag = .true.
       case Default
          call MessageNotify('E', module_name, &
               & "The specified type of budget analysis '%c' is invalid.", c1=trim(budgetAnaName(n)) )
       end select
    end do

    call prepair_Output(diagVar_gthsInfo)

  end subroutine BudgetAnalysis_Init

  !>
  !!
  !!
  subroutine BudgetAnalysis_Final()

    ! 実行文; Executable statements
    !

    if( energyBudgAnaFlag ) then
       call HistoryClose(hst_energyBudget)
    end if

    if( angMomBudgAnaFlag ) then
       call HistoryClose(hst_angMomBudget)
    end if

  end subroutine BudgetAnalysis_Final

  subroutine BudgetAnalysis_perform( &
       & xyz_U, xyz_V, xyz_SigDot, xyz_DensEdd, &
       & xy_SurfPress, xy_totDepth )

    ! 宣言文; Declare statements
    !
    real(DP), intent(in), dimension(0:iMax-1, jMax, 0:kMax) :: &
         & xyz_U, xyz_V, xyz_SigDot, xyz_DensEdd
    real(DP), intent(in), dimension(0:iMax-1, jMax) :: xy_SurfPress, xy_totDepth
 
    ! 実行文; Executable statements
    !
    
    if( energyBudgAnaFlag ) call analyze_energyBudget()
    if( angMomBudgAnaFlag ) call analyze_angMomBudget()

  contains
    subroutine analyze_energyBudget()

use wa_module
use HydroBouEqSolverRHS_mod

      real(DP) :: KEnAvg, PEnAvg
      
      real(DP) :: PEConvert, HDiffDisp, HDiffDisp2, VDiffDisp, VDiffDisp2
      real(DP) :: KEInput, KEInput2, NoHDivCorrecTmp, NoHDivCorrec, NoHDivCorrec2
      real(DP) :: Adv
      real(DP) :: xyz_DuDSig(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_DvDSig(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_DPsiDSig(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_DChiDSig(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_CosLat(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_Div(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_Vor(0:iMax-1,jMax,0:kMax)
      real(DP) :: w_Zero(lMax,0:kMax)

 real(DP) :: xy(0:iMax-1,jMax), xyz_DwDz(0:iMax-1,jMax,0:kMax)
 real(DP) :: wz_DivRHS(lMax,0:kMax)
 real(DP) :: wz_VorRHS(lMax,0:kMax)
 real(DP) :: xyz_AdvPsi(0:iMax-1,jMax,0:kMax), xyz_AdvChi(0:iMax-1,jMax,0:kMax)
 real(DP), parameter :: dt = 36000d0

      call MessageNotify("M", module_name, "Analyze energy budget..")
      
      !
      KEnAvg = eval_kineticEnergyAvg(xyz_U, xyz_V)
      PEnAvg = eval_potentialEnergyAvg(xyz_DensEdd, xy_totDepth) 
      call HistoryPut( BUDGETANAKEY_KEAVG, KEnAvg, hst_energyBudget )
      call HistoryPut( BUDGETANAKEY_PEAVG, PEnAvg, hst_energyBudget )
      call HistoryPut( BUDGETANAKEY_TEAVG, KEnAvg + PEnAvg, hst_energyBudget )

      !
      xyz_Div = eval_Div(xyz_U, xyz_V)
      xyz_Vor = eval_Vor(xyz_U, xyz_V)

      xyz_DuDSig = xyz_wt( wt_DSig_wt(wt_xyz(xyz_U)) ) 
      xyz_DvDSig = xyz_wt( wt_DSig_wt(wt_xyz(xyz_V)) ) 
      xyz_DPsiDSig = xyz_wt( wa_LaplaInv_wa( wt_DSig_wt(wt_xyz(xyz_Vor)) )*RPlanet**2 ) 
      xyz_DChiDSig = xyz_wt( wa_LaplaInv_wa( wt_DSig_wt(wt_xyz(xyz_Div)) )*RPlanet**2 ) 

      xyz_CosLat = cos(xyz_Lat)

      PEConvert = -AvrLonLat_xy( xy_totDepth * xy_IntSig_BtmToTop_xyz( (xyz_DensEdd/RefDens)*Grav*xyz_SigDot ) )

      KEInput2 = vDiffCoef*AvrLonLat_xy( &
           & (xyz_U(:,:,0)*xyz_DuDSig(:,:,0) + xyz_V(:,:,0)*xyz_DvDSig(:,:,0))/xy_totDepth**2 &
           & )
      VDiffDisp2 = - vDiffCoef* AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
           & xyz_DuDSig**2 + xyz_DvDSig**2 &
           & )/xy_totDepth**2 )

      xyz_DuDSig = xyz_CosLat*xyz_AlphaOptr_wz(wz_xyz(xyz_DChiDSig), -wz_xyz(xyz_DPsiDSig))
      xyz_DvDSig = xyz_CosLat*xyz_AlphaOptr_wz(wz_xyz(xyz_DPsiDSig),  wz_xyz(xyz_DChiDSig))
      KEInput = vDiffCoef*AvrLonLat_xy(  &
           &   (xyz_U(:,:,0) * xyz_DuDSig(:,:,0) + xyz_V(:,:,0) * xyz_DvDSig(:,:,0))/xy_totDepth**2 & 
           & )
      
      VDiffDisp = vDiffCoef* AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
           &    xyz_U*xyz_wt(wt_DSig_wt(wt_xyz(xyz_DuDSig))) &
           &  + xyz_V*xyz_wt(wt_DSig_wt(wt_xyz(xyz_DvDSig))) &
           & )/xy_totDepth**2 ) - KEInput

      HDiffDisp2 = hDiffCoef*AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
           &    xyz_U * ( xyz_wz(wz_Lapla2D_wz(wz_xyz(xyz_U))) - xyz_U/(RPlanet*xyz_CosLat)**2 ) &
           &  + xyz_V * ( xyz_wz(wz_Lapla2D_wz(wz_xyz(xyz_V))) - xyz_V/(RPlanet*xyz_CosLat)**2 ) &
           & ))
     
      HDiffDisp = hDiffCoef*AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
           &   xyz_U * xyz_CosLat*xyz_AlphaOptr_wz(wz_xyz(xyz_Div), -wz_xyz(xyz_Vor)) &
           & + xyz_V * xyz_CosLat*xyz_AlphaOptr_wz(wz_xyz(xyz_Vor),  wz_xyz(xyz_Div)) &
           & ) )

      w_Zero = 0d0      
      NoHDivCorrec2 = AvrLonLat_xy(  &
           & xy_IntSig_BtmToTop_xyz(xyz_Div)*xy_SurfPress/RefDens &
           & )
      NoHDivCorrec = - AvrLonLat_xy(  &
           &   xy_IntSig_BtmToTop_xyz(xyz_U*xyz_CosLat)*xy_AlphaOptr_w(w_xy(xy_SurfPress/RefDens), w_Zero) &
           & + xy_IntSig_BtmToTop_xyz(xyz_V*xyz_CosLat)*xy_AlphaOptr_w(w_Zero, w_xy(xy_SurfPress/RefDens)) &
           &  )

!!

      call calc_VorEqDivEqInvisRHS(wz_VorRHS, wz_DivRHS, &
         & xyz_Vor, xyz_U*xyz_CosLat, xyz_V*xyz_CosLat, xy_w(w_Zero), xyz_DensEdd, &
         & Diagnose_PressBaroc(xy_totDepth, xyz_DensEdd), Diagnose_GeoPot(xy_totDepth), xyz_SigDot)

      Adv = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
           &    xyz_U * RPlanet**2 * xyz_CosLat*xyz_AlphaOptr_wz(wa_LaplaInv_wa(wz_DivRHS), -wa_LaplaInv_wa(wz_VorRHS)) &
           &  + xyz_V * RPlanet**2 * xyz_CosLat*xyz_AlphaOptr_wz(wa_LaplaInv_wa(wz_VorRHS),  wa_LaplaInv_wa(wz_DivRHS))  &
           & ) )


!!
      call HistoryPut( BUDGETANAKEY_PE2KEAVG, PEConvert, hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_KE2PEAVG, -PEConvert, hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_KEINPUTSURFAVG, KEInput, hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_HDIFFDISPAVG, HDiffDisp, hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_VDIFFDISPAVG, VDiffDisp, hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_SURFPRESSWORKAVG, NoHDivCorrec, hst_energyBudget)
      call HistoryPut( BUDGETANAKEY_KEGENNETAVG, & 
           & PEConvert + KEInput + VDiffDisp + HDiffDisp + NoHDivCorrec, &
           & hst_energyBudget)

!!$write(*,*) RefDens*vDiffCoef*xyz_DuDSig(1,:,0)/xy_totDepth(1,1)a
xy = xy_IntSig_BtmToTop_xyz(xyz_Div)
xyz_DwDz = xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot)))
write(*,*) "Check Continuity:"
write(*,*) (xyz_DwDz(1,1:16,0)+xyz_Div(1,1:16,0))/maxval(xyz_Div)
write(*,*) "="
write(*,*) "HDiffDisp:", HDiffDisp2, HDiffDisp
write(*,*) "VDiffDisp:", VDiffDisp2, VDiffDisp
write(*,*) "NoHDivCorrec:", NoHDivCorrec2, NoHDivCorrec
write(*,*) "KEInput:", KEInput2, KEInput
write(*,*) "Adv:", Adv

    end subroutine analyze_energyBudget

    subroutine analyze_angMomBudget()

      real(DP) :: AngMomAvg

      call MessageNotify("M", module_name, "Analyze angular momentum budget..")

      AngMomAvg = eval_angularMomAvg(xyz_U)
      call HistoryPut( BUDGETANAKEY_ANGMOMAVG, AngMomAvg, hst_angMomBudget )

    end subroutine analyze_angMomBudget

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

    if ( energyBudgAnaFlag ) then
       call HistoryCreate( & 
            & file= trim(diagVar_gthsInfo%FilePrefix) // 'EnergyBudget.nc', title='energy budget analysis', &
            & source='Dennou-OGCM', &
            & institution='Dennou-OGCM project', &
            & dims=(/'t'/), dimsizes=(/ 0 /), &
            & longnames=(/'time'/),&
            & units=(/ diagVar_gthsInfo%intUnit /), &
            & origin=real(0), interval=real(diagVar_gthsInfo%intValue), &
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
             & varname=BUDGETANAKEY_KEGENNETAVG, dims=(/'t'/), &
             & longname='global mean of net kinetic energy generation', &
             & units='J*m-3*kg-1*s-1', xtype='double',&
             & history=hst_energyBudget)

        call HistoryAddVariable( & 
             & varname=BUDGETANAKEY_SURFPRESSWORKAVG, dims=(/'t'/), &
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

