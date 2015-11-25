!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SGSConvAdjust_mod 

  ! モジュール引用; Use statements
  !
  !
  use dc_types, only: DP, TOKEN
  use dc_message, only: MessageNotify

  use Constants_mod, only: &
       & RPlanet

  use SpmlUtil_mod

  use DiagnoseUtil_mod

  use EOSDriver_mod

  use GridSet_mod, only: &
       & iMax, jMax, kMax, &
       & z_LyrThickSig

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  interface SGSConvAdjust_perform
     module procedure SGSConvAdjust_perform_1D
     module procedure SGSConvAdjust_perform_GCMDriver
  end interface SGSConvAdjust_perform

  public :: SGSConvAdjust_Init, SGSConvAdjust_Final
  public :: SGSConvAdjust_perform
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SGSConvAdjust_mod' !< Module Name

  real(DP) :: STABILITY_THRESHOLD 
  
contains

  !>
  !!
  !!
  subroutine SGSConvAdjust_Init()

    ! 実行文; Executable statements
    !

    STABILITY_THRESHOLD = 0d0

    
  end subroutine SGSConvAdjust_Init

  !>
  !!
  !!
  subroutine SGSConvAdjust_Final()

    ! 実行文; Executable statements
    !

  end subroutine SGSConvAdjust_Final

  !> @brief 
  !!
  !!
  subroutine SGSConvAdjust_perform_GCMDriver( &
       & xyz_PTemp, xyz_Salt, xy_totDepth, &
       & xyz_isAdjustOccur )

    !
    use SpmlUtil_mod, only: g_Sig

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: &
         & xyz_PTemp, xyz_Salt
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    logical, intent(inout) :: xyz_isAdjustOccur(0:iMax-1,jMax, 0:kMax)

    ! 局所変数
    ! Local variables
    !
    integer :: i, j
    real(DP), dimension(0:kMax) ::  z_PTemp, z_Salt
    
    ! 実行文; Executable statement
    !


    !$omp parallel do private(i, z_PTemp, z_Salt)
    do j=1,jMax
       do i=0,iMax-1
          z_PTemp(:) = xyz_PTemp(i,j,:)
          z_Salt(:) = xyz_Salt(i,j,:)

          call SGSConvAdjust_perform(z_PTemp, z_Salt,      & !(inout)
               & xy_totDepth(i,j)*z_LyrThickSig(:),        & !(in)
               & xyz_isAdjustOccur(i,j,:) )

          xyz_PTemp(i,j,:) = z_PTemp(:)
          xyz_Salt(i,j,:) = z_Salt(:)
       end do
    end do
    
  end subroutine SGSConvAdjust_perform_GCMDriver

  subroutine SGSConvAdjust_perform_1D( z_PTemp, z_Salt, &
       & z_LyrThick, isAdjustOccur )

    use Constants_mod
    
    ! 実行文; Executable statement
    ! 
    real(DP), dimension(0:kMax), intent(inout) :: z_PTemp, z_Salt
    real(DP), dimension(0:kMax), intent(in) :: z_LyrThick
    logical, dimension(0:kMax), intent(inout) :: isAdjustOccur

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:kMax) :: z_DensPotEdd, z_PressRef
    integer :: k_, k
    logical :: isColumnStable
    integer :: nItr

    integer :: r
    real(DP) :: r_RefPress(1:kMax)
    logical :: r_StaticStableFlag(1:kMax)
    real(DP) :: DensPotPair(2)
    
    ! 実行文; Executable statement
    !

    !
    z_PressRef = 0d0

    !
    call EOSDriver_Eval(rhoEdd=z_DensPotEdd,        & ! (out)
         & theta=z_PTemp, S=z_Salt, p=z_PressRef    & ! (in)
         & )

    r_StaticStableFlag(:) = .true.
    do k_=1, kMax
       r_RefPress(k_) = sum(z_LyrThick(0:k_-1))*RefDens*Grav
       call EOSDriver_Eval(rhoEdd=DensPotPair(1), &
            & theta=z_PTemp(k_-1), S=z_Salt(k_-1), p=r_RefPress(k_))
       call EOSDriver_Eval(rhoEdd=DensPotPair(2), &
            & theta=z_PTemp(k_), S=z_Salt(k_), p=r_RefPress(k_))

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
       !
!!$       do k_=0,kMax-1
!!$          if(z_DensPotEdd(k_) > z_DensPotEdd(k_+1)) then
!!$             call mix_below1stUnstableLyr(z_DensPotEdd, k_)
!!$             exit
!!$          end if
!!$          if(k_==kMax-1) then
!!$             ! If k has reached kMax-1, the sea water column is stable, and exit this subroutine.
!!$             isColumnStable = .true.
!!$          end if
!!$       end do
!!$
       do k_=1, kMax
          if(.not. r_StaticStableFlag(k_)) then
             call mix_UnstableLyr(k_)
             exit
          end if
          if(k_==kMax) isColumnStable = .true.
       end do
    end do

    !
    !
  contains
    subroutine mix_UnstableLyr(k_)
      integer, intent(in) :: k_
      integer :: nLyrMix
      real(DP) :: DensPotPair(2)
      integer :: LyrLId, LyrUId, m

      
      nLyrMix = 1
      do m=2, kMax-k_+2
         LyrLId = k_-1
         LyrUId = k_+m-2
         nLyrMix = nLyrMix + 1

         call vmixing_PTempAndSalt(z_PTemp(LyrLId:LyrUId), z_Salt(LyrLId:LyrUId), &
              & z_LyrThick(LyrLId:LyrUId), nLyrMix )

         if(LyrUId+1 /= kMax+1) then
            call EOSDriver_Eval(rhoEdd=DensPotPair(:), &
                 & theta=z_PTemp(LyrUId:LyrUId+1), S=z_Salt(LyrUId:LyrUId+1), p=spread(r_RefPress(k_+nLyrMix-1),1,2) )
            if(DensPotPair(1) <= DensPotPair(2)) exit
         end if
      end do

      r_StaticStableFlag(k_:k_+nLyrMix-2) = .true.
      
      if(LyrLId /= 0) then
         call EOSDriver_Eval(rhoEdd=DensPotPair(:), &
              & theta=z_PTemp(LyrLId-1:LyrLId), S=z_Salt(LyrLId-1:LyrLId), p=spread(r_RefPress(k_-1),1,2) )
            if(DensPotPair(1) > DensPotPair(2)) r_StaticStableFlag(k_-1) = .false.
      end if

      isAdjustOccur(LyrLId:LyrUId) = .true.

    end subroutine mix_UnstableLyr
    
    subroutine mix_below1stUnstableLyr(z_DensPotEdd, k)

      real(DP), intent(inout) :: z_DensPotEdd(0:kMax)
      integer, intent(in) :: k

      integer :: nLyrMix

      isAdjustOccur(k) = .true.

      do nLyrMix=2,kMax-k+1
         isAdjustOccur(k+nLyrMix-1) = .true.
         call vmixing_PTempAndSalt(z_PTemp(k:k+nLyrMix-1), z_Salt(k:k+nLyrMix-1), &
              & z_LyrThick(k:k+nLyrMix-1), nLyrMix )

         call EOSDriver_Eval( rhoEdd=z_DensPotEdd(k:k+nLyrMix-1), &
              & theta=z_PTemp(k:k+nLyrMix-1), S=z_Salt(k:k+nLyrMix-1), p=z_PressRef(k:k+nLyrMix-1) &
              & )

         if(k+nLyrMix==kMax+1) exit
         if(z_DensPotEdd(k+nLyrMix-1) <= z_DensPotEdd(k+nLyrMix)) exit
      end do

    end subroutine mix_below1stUnstableLyr

  end subroutine SGSConvAdjust_perform_1D

  subroutine vmixing_PTempAndSalt(PTemp, Salt, LyrThick, nLayerMix)
    integer, intent(in) :: nLayerMix
    real(DP), intent(inout) :: PTemp(nLayerMix)
    real(DP), intent(inout) :: Salt(nLayerMix)
    real(DP), intent(in) :: LyrThick(nLayerMix)

    real(DP) :: mixed_PTemp, mixed_Salt
    real(DP) :: totLyrThick

    totLyrThick = sum(LyrThick(:))
    mixed_PTemp = sum(PTemp*LyrThick(:))/totLyrThick
    mixed_Salt = sum(Salt*LyrThick(:))/totLyrThick

    PTemp(:) = mixed_PTemp
    Salt(:) = mixed_Salt

  end subroutine vmixing_PTempAndSalt

end module SGSConvAdjust_mod

