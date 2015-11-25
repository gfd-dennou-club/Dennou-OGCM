!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SGSSlowConvAdjust_mod 

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

  use GridSet_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  interface SGSSlowConvAdjust_perform
     module procedure SGSSlowConvAdjust_perform_1D
     module procedure SGSSlowConvAdjust_perform_GCMDriver
  end interface SGSSlowConvAdjust_perform

  public :: SGSSlowConvAdjust_Init, SGSSlowConvAdjust_Final
  public :: SGSSlowConvAdjust_perform
  
  ! 非公開手続き
  ! Private procedure
  !

  real(DP), public, parameter :: DEFAULT_MIXTIME = 12d0*3600d0
  real(DP), public, parameter :: DEFAULT_CONVLAYER_DEPTH = 1d3

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SGSSlowConvAdjust_mod' !< Module Name

  real(DP) :: DelTime, MixTime, ConvLyrDepth
  real(DP) :: STABILITY_THRESHOLD 
  real(DP) :: hMixForDelTime


contains

  !>
  !!
  !!
  subroutine SGSSlowConvAdjust_Init(GCMTimeStep, ConvMixTime, ConvLayerDepth)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: GCMTimeStep
    real(DP), intent(in), optional :: ConvMixTime, ConvLayerDepth

    ! 実行文; Executable statements
    !
    
    DelTime = GCMTimeStep
    
    !
    MixTime = DEFAULT_MIXTIME
    ConvLyrDepth = DEFAULT_CONVLAYER_DEPTH
    if(present(ConvMixTime)) MixTime = ConvMixTime
    if(present(ConvLayerDepth)) ConvLyrDepth = ConvLayerDepth

    !
    hMixForDelTime = 0.4d0*ConvLyrDepth*sqrt(DelTime/MixTime)
    
    call MessageNotify('M', module_name, &
         & "Set mixing length during one time step in OGCM.. %f [m]", d=(/ hMixForDelTime /))

  end subroutine SGSSlowConvAdjust_Init

  !>
  !!
  !!
  subroutine SGSSlowConvAdjust_Final()

    ! 実行文; Executable statements
    !

  end subroutine SGSSlowConvAdjust_Final

  !> @brief 
  !!
  !!
  subroutine SGSSlowConvAdjust_perform_GCMDriver( &
       & xyz_PTemp, xyz_Salt, xy_totDepth, &
       & xyz_isAdjustOccur )

    !
    use SpmlUtil_mod, only: g_Sig

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: &
         & xyz_PTemp, xyz_Salt
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    logical, intent(out) :: xyz_isAdjustOccur(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    integer :: i, j
    real(DP), dimension(0:kMax) ::  z_PTemp, z_Salt, z_Depth
    
    ! 実行文; Executable statement
    !


    !$omp parallel do private(i, z_PTemp, z_Salt, z_Depth)
    do j=1,jMax
       do i=0,iMax-1

          z_PTemp(:) = xyz_PTemp(i,j,:)
          z_Salt(:) = xyz_Salt(i,j,:)
          z_Depth(:) = xy_totDepth(i,j)*g_Sig(:)

          call SGSSlowConvAdjust_perform(z_PTemp, z_Salt, & !(inout)
               & z_Depth,                             & !(in)
               & xyz_isAdjustOccur(i,j,:) )

          xyz_PTemp(i,j,:) = z_PTemp(:)
          xyz_Salt(i,j,:) = z_Salt(:)

       end do
    end do
    

  end subroutine SGSSlowConvAdjust_perform_GCMDriver

  subroutine SGSSlowConvAdjust_perform_1D( z_PTemp, z_Salt, &
       & z_Depth, isAdjustOccur )

    ! 実行文; Executable statement
    ! 
    real(DP), dimension(0:kMax), intent(inout) :: z_PTemp, z_Salt
    real(DP), dimension(0:kMax), intent(out) :: z_Depth
    logical, intent(out) :: isAdjustOccur(0:kMax)

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:kMax) :: z_DensPotEdd, z_PressRef, z_LyrThick

    real(DP) :: A(kMax,0:kMax), B(kMax+1), Work(2*(kMax+1))
    integer :: info
    
    integer :: nLayer
    integer :: k_, k, n_, n, N_homo
    integer :: z_NMix(0:kMax)


    ! 実行文; Executable statement
    !

    !
    z_PressRef = 0d0

    !
    !
    A = 0d0; B = 0d0
    do k=1,kMax
       A(k,k-1:k) = 0.5d0
       B(k) = z_Depth(k-1) - z_Depth(k) !g_Sig(k-1) - g_Sig(k)
    end do
    A(1,0) = 1; A(kMax,kMax) = 1

    call DGELS('N', kMax, kMax+1, 1, A, kMax, B, kMax+1, Work, 2*(kMax+1), info)
    z_LyrThick = b


    !
    !
    do k=0, kMax-1
       do n_=1,kMax-k
          if(sum(z_LyrThick(k+1:k+n_)) >= hMixForDelTime .or. k+n_==kMax) then
             n = n_; exit
          end if
       end do

       if( abs(sum(z_LyrThick(k+1:k+n))-hMixForDelTime) > abs(sum(z_LyrThick(k+1:k+n-1))-hMixForDelTime)) then
          n = n - 1
       end if
       z_NMix(k) = n
!       write(*,*) k, sum(z_LyrThick(k+1:k+n)), z_Depth(k), n
    end do
    z_NMix(kMax) = 0
    
    !
    call EOSDriver_Eval(rhoEdd=z_DensPotEdd,        & ! (out)
         & theta=z_PTemp, S=z_Salt, p=z_PressRef    & ! (in)
         & )


    !
    isAdjustOccur(:) = .false.
    
    !
    k_ = kMax-1
    do while(k_ >= 0)
       if(z_DensPotEdd(k_) > z_DensPotEdd(k_+1)) then
          isAdjustOccur(k_:k_+1) = .true.
          call mix_below1stUnstableLyr(z_DensPotEdd, k_, N_homo)

          k_ = k_ - N_homo
       else
          k_ = k_ - 1
       end if
    end do


    !
    !
  contains
    subroutine mix_below1stUnstableLyr(z_DensPotEdd, k, N_h)
      real(DP), intent(in) :: z_DensPotEdd(0:kMax)
      integer, intent(in) :: k
      integer, intent(out) :: N_h

      real(DP), parameter :: EPS = 1d-16
      integer :: N_Mix, k_
      
      N_h = 1
      if(k /= 0) then
         do while( abs(z_DensPotEdd(k-N_h)-z_DensPotEdd(k)) < EPS )
            N_h = N_h + 1
            if(k-N_h-1 <= 0) exit
         end do
      end if

      N_Mix = 1
      do k_=k, k+z_NMix(k)
         if(z_DensPotEdd(k_) < z_DensPotEdd(k_+1) .or. k+N_Mix==kMax) exit
         N_Mix = N_Mix + 1
      end do


      call vmixing_Tracer(z_PTemp, z_LyrThick, k, N_h, N_Mix)
      call vmixing_Tracer(z_Salt, z_LyrThick, k, N_h, N_Mix)


    end subroutine mix_below1stUnstableLyr

  end subroutine SGSSlowConvAdjust_perform_1D

  subroutine vmixing_Tracer(Tracer, LyrThick, k, N_h, N_mix)
    integer, intent(in) :: k, N_h, N_mix
    real(DP), intent(inout) :: Tracer(0:kMax)
    real(DP), intent(in) :: LyrThick(0:kMax)

    real(DP) :: delTracer
    real(DP) :: w
    integer :: k_u, k_l
        
    k_u = k - N_h + 1
    k_l = k + N_mix
   
    w = sum(LyrThick(k_u:k))/sum(LyrThick(k_u:k_l))

!!$if(k_u < 0 .or. k_l > kMax) then
!!$   write(*,*) "k=", k, ", range=", k_u, k_l, ", nmix=", N_mix
!!$stop
!!$end if
    delTracer = Tracer(k) - Tracer(k+1)
    Tracer(k_u:k) = Tracer(k_u:k) - (1d0 - w)*delTracer
    Tracer(k+1:k_l) = Tracer(k+1:k_l) + w*delTracer

  end subroutine vmixing_Tracer

end module SGSSlowConvAdjust_mod

