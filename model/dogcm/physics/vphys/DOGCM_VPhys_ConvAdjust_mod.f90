!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_VPhys_ConvAdjust_mod 

  ! モジュール引用; Use statements
  !
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM
  
  use DOGCM_Admin_Constants_mod, only: &
       & RPlanet, RefDens, Grav

  use EOSDriver_mod, only: &
       & EOSDriver_Eval

  use DOGCM_Admin_Grid_mod, only: &
       & iMax, jMax,              &
       & KS, KE

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  interface DOGCM_VPhys_ConvAdjust_perform
     module procedure DOGCM_VPhys_ConvAdjust_perform_1D
  end interface DOGCM_VPhys_ConvAdjust_perform

  public :: DOGCM_VPhys_ConvAdjust_Init, DOGCM_VPhys_ConvAdjust_Final
  public :: DOGCM_VPhys_ConvAdjust_AddMixingTerm

  public :: DOGCM_VPhys_ConvAdjust_perform


  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_VPhys_ConvAdjust_mod' !< Module Name

  real(DP) :: STABILITY_THRESHOLD 

  integer, parameter :: ExceptBC_OFS = 0
  
contains

  !>
  !!
  !!
  subroutine DOGCM_VPhys_ConvAdjust_Init()

    ! 実行文; Executable statements
    !

    STABILITY_THRESHOLD = 0d0

    
  end subroutine DOGCM_VPhys_ConvAdjust_Init

  !>
  !!
  !!
  subroutine DOGCM_VPhys_ConvAdjust_Final()

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_VPhys_ConvAdjust_Final

  !> @brief 
  !!
  !!
  subroutine DOGCM_VPhys_ConvAdjust_AddMixingTerm( &
       & xyz_PTemp_RHS, xyz_Salt_RHS,                     & ! (inout)
       & xyz_ConvIndex,                                   & ! (inout)       
       & xyz_PTemp, xyz_Salt,                             & ! (in)
       & xyz_Z, z_LyrThickIntWt, dt                       & ! (in)
       & )

    use SpmlUtil_mod

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xyz_PTemp_RHS(0:iMax-1,jMax,KS:KE)
    real(DP), intent(inout) :: xyz_Salt_RHS(0:iMax-1,jMax,KS:KE)
    real(DP), intent(inout) :: xyz_ConvIndex(0:iMax-1,jMax,KS:KE)    
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: z_LyrThickIntWt(KS:KE)
    real(DP), intent(in) :: dt
    
    ! 局所変数
    ! Local variables
    !
    integer :: i, j
    real(DP) :: z_PTemp(KS:KE)
    real(DP) :: z_PTempOri(KS:KE)
    real(DP) :: z_Salt(KS:KE)
    real(DP) :: z_Z(KS:KE)
    logical  :: z_isAdjustOccur(KS:KE)

    real(DP) :: PTempVInt

    ! 実行文; Executable statement
    !

    
    !$omp parallel do private(i, PTempVInt, z_PTempOri, z_PTemp, z_Salt, z_Z, z_isAdjustOccur)
    do j=1,jMax
       do i=0,iMax-1
          z_PTemp(:) = xyz_PTemp(i,j,:)
          z_Salt(:) = xyz_Salt(i,j,:)
          z_Z(:) = xyz_Z(i,j,:)

          z_PTempOri = z_PTemp
          PTempVInt = sum(z_LyrThickIntWt(:)*z_PTemp)
          
          call DOGCM_VPhys_ConvAdjust_perform_1D(z_PTemp, z_Salt,      & ! (inout)
               & z_Z, z_LyrThickIntWt,                        & ! (in)
               & z_isAdjustOccur    )                           ! (out)
          
          xyz_PTemp_RHS(i,j,:) = xyz_PTemp_RHS(i,j,:) + (z_PTemp(:) - xyz_PTemp(i,j,:))/dt
          xyz_Salt_RHS(i,j,:) = xyz_Salt_RHS(i,j,:) + (z_Salt(:) - xyz_Salt(i,j,:))/dt

          if( (PTempVInt - sum(z_LyrThickIntWt(:)*z_PTemp))/PTempVInt > 1d-12 ) then
             write(*,*) "before=", PTempVInt, "after=", sum(z_LyrThickIntWt(:)*z_PTemp)
             write(*,*) "LyrThickWt=", z_LyrThickIntWt(:)
             write(*,*) "z_PTemp=", z_PTemp
             call MessageNotify('E', module_name, &
                  & "Integrated potential temperature over a column is not conserved. Check!") 
          end if
          

          where( z_isAdjustOccur(:) )
             xyz_ConvIndex(i,j,:) = xyz_ConvIndex(i,j,:) + 1d0
          end where
          
       end do
    end do
    
  end subroutine DOGCM_VPhys_ConvAdjust_AddMixingTerm
  
  !----------------------------------------------------------
  
  subroutine DOGCM_VPhys_ConvAdjust_perform_1D( z_PTemp, z_Salt, & ! (inout) 
       & z_Z, z_LyrThickIntWt, isAdjustOccur )                     ! (in)
    
    ! 実行文; Executable statement
    ! 
    real(DP), intent(inout) :: z_PTemp(KS:KE)
    real(DP), intent(inout) :: z_Salt(KS:KE)
    real(DP), intent(in) :: z_Z(KS:KE)
    real(DP), intent(in) :: z_LyrThickIntWt(KS:KE)
    logical, intent(out) :: isAdjustOccur(KS:KE)

    ! 局所変数
    ! Local variables
    !
    integer :: k_
    integer :: k
    logical :: isColumnStable
    integer :: nItr

    integer :: r
    real(DP) :: r_RefPress(KS-1:KE)
    logical :: r_StaticStableFlag(KS-1:KE)
    real(DP) :: DensPotPair(2)
    
    ! 実行文; Executable statement
    !

    !
    !
    r_StaticStableFlag(:) = .true.
    do k_=KS, KE-1
       r_RefPress(k_) = - RefDens*Grav*0.5d0*(z_Z(k_) + z_Z(k_+1))

       call EOSDriver_Eval( rhoEdd=DensPotPair(:),       & ! (out)
            & theta=z_PTemp(k_:k_+1), S=z_Salt(k_:k_+1), & ! (in)
            & p=(/ r_RefPress(k_), r_RefPress(k_) /) )     ! (in)

       if(DensPotPair(1) > DensPotPair(2)) then
          r_StaticStableFlag(k_) = .false.
       end if
    end do
    
    !
    k_ = 0
    nItr = 0
    isColumnStable = .false.
    isAdjustOccur(:) = .false.

    
    !
    do while(.not. isColumnStable)
       nItr = nItr + 1
!write(*,*) "itr=", nItr, "*densPot:", z_DensPotEdd

       do k_ = KS, KE-1

          if(.not. r_StaticStableFlag(k_)) then
             call mix_UnstableLyr(k_)
             exit
          end if

          if(k_ == KE-1) then
             isColumnStable = .true.
          end if
          
       end do
       
    end do

    !
    !
  contains
    subroutine mix_UnstableLyr(r_)
      integer, intent(in) :: r_
      integer :: nMixPair
      real(DP) :: DensPotPair(2)
      integer :: LyrLId
      integer :: LyrUId
      integer :: m

      
      nMixPair = 0
      LyrLId = r_
      
      do m=2, KE-r_+2

         LyrUId = r_ + m - 2
         nMixPair = nMixPair + 1

         call vmixing_PTempAndSalt(z_PTemp(LyrLId:LyrUId), z_Salt(LyrLId:LyrUId), & ! (inout)
              & z_LyrThickIntWt(LyrLId:LyrUId) )                                    ! (in)
         
         if(LyrUId /= KE) then
            call EOSDriver_Eval(rhoEdd=DensPotPair(:),                        &   ! (out)
                 & theta=z_PTemp(LyrUId:LyrUId+1), S=z_Salt(LyrUId:LyrUId+1), &   ! (in)
                 & p=(/ r_RefPress(r_+nMixPair), r_RefPress(r_+nMixPair) /) )     ! (in)
            
            if (DensPotPair(1) < DensPotPair(2)) exit
         end if
      end do
      r_StaticStableFlag(r_:r_+nMixPair-1) = .true.
      
      if(LyrLId /= KS) then
         call EOSDriver_Eval(rhoEdd=DensPotPair(:),                        & ! (out)
              & theta=z_PTemp(LyrLId-1:LyrLId), S=z_Salt(LyrLId-1:LyrLId), & ! (in)
              & p=(/ r_RefPress(r_-1), r_RefPress(r_-1) /) )                 ! (in)
         
            if(DensPotPair(1) > DensPotPair(2)) r_StaticStableFlag(r_-1) = .false.
      end if

      !
      isAdjustOccur(LyrLId:LyrUId) = .true.

    end subroutine mix_UnstableLyr
    
!!$    subroutine mix_below1stUnstableLyr(z_DensPotEdd, k)
!!$
!!$      real(DP), intent(inout) :: z_DensPotEdd(0:kMax)
!!$      integer, intent(in) :: k
!!$
!!$      integer :: nLyrMix
!!$
!!$      isAdjustOccur(k) = .true.
!!$
!!$      do nLyrMix=2,kMax-k+1
!!$         isAdjustOccur(k+nLyrMix-1) = .true.
!!$         call vmixing_PTempAndSalt(z_PTemp(k:k+nLyrMix-1), z_Salt(k:k+nLyrMix-1), &  ! (inout)
!!$              & z_LyrThickIntWt(k:k+nLyrMix-1), nLyrMix                           &  ! (in)
!!$              & )
!!$
!!$         call EOSDriver_Eval( rhoEdd=z_DensPotEdd(k:k+nLyrMix-1),                 &  ! (out)
!!$              & theta=z_PTemp(k:k+nLyrMix-1), S=z_Salt(k:k+nLyrMix-1),            &  ! (in)
!!$              & p=z_PressRef(k:k+nLyrMix-1) )                                        ! (in)
!!$
!!$         if(k+nLyrMix==kMax+1) exit
!!$         if(z_DensPotEdd(k+nLyrMix-1) <= z_DensPotEdd(k+nLyrMix)) exit
!!$      end do
!!$
!!$    end subroutine mix_below1stUnstableLyr

  end subroutine DOGCM_VPhys_ConvAdjust_perform_1D

  subroutine vmixing_PTempAndSalt(PTemp, Salt, LyrThick)
    real(DP), intent(in) :: LyrThick(:)
    real(DP), intent(inout) :: PTemp(size(LyrThick))
    real(DP), intent(inout) :: Salt(size(LyrThick))

    real(DP) :: mixed_PTemp
    real(DP) :: mixed_Salt
    real(DP) :: totLyrThick
    
    totLyrThick = sum(LyrThick(:))
    mixed_PTemp = sum(PTemp*LyrThick(:))/totLyrThick
    mixed_Salt = sum(Salt*LyrThick(:))/totLyrThick

    PTemp(:) = mixed_PTemp
    Salt(:) = mixed_Salt

  end subroutine vmixing_PTempAndSalt

end module DOGCM_VPhys_ConvAdjust_mod

