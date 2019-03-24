program SeaWaterProp_mod_test
  use dc_types
  use dc_test

  implicit none

  call test_FM83()

contains
  subroutine test_FM83()
    
    ! モジュール引用; Use statement
    !
    use SeaWaterProp_FM83_mod, only: &
         & SeaWaterProp_FM83_Init, SeaWaterProp_FM83_Final, &
         & SeaWaterProp_FM83_AdLapseRate, SeaWaterProp_FM83_Temp2PTemp, &
         & SeaWaterProp_FM83_HeatCapacity

    ! 作業変数
    ! Work variable
    !
    real(DP), parameter :: answer_AdLapseRate = 3.255976d-04
    real(DP), parameter :: answer_InSituTemp2Ptemp1 =  36.89073d0
    real(DP), parameter :: answer_InSituTemp2Ptemp2 = - 1.09749d0
    real(DP), parameter :: ansert_HeatCapacity = 3849.500d0
 
    real(DP), parameter :: acceptableError = 1d-05   ! The error relative to true value is less than 0.001%. 
real(DP) :: cp
    ! 実行文; Executable statement
    !

    call SeaWaterProp_FM83_Init()

    call  AssertLessThan( &
         &  message="Check SeaWaterProp_FM83_AdLapseRate", &
         &  answer=acceptableError, check= abs( answer_AdLapseRate - &
         &     SeaWaterProp_FM83_AdLapseRate( &
         &       p   = 1d04, &   ! (in)  [dbar]
         &       S   = 40d0, &   ! (in)  [psu]
         &       t   = 40d0 &    ! (in)  [degree C]
         &      ) &              ! (out) [degree C]
         &   ) / answer_AdLapseRate &
         & )
    
     call  AssertLessThan( &
         &  message="Check SeaWaterProp_FM83_Temp2PTemp", &
         &  answer=acceptableError, check= abs( answer_InSituTemp2Ptemp1 - &
         &     SeaWaterProp_FM83_Temp2PTemp( &
         &       P0   = 1d04, &  ! (in)  [dbar]
         &       S0   = 40d0, &  ! (in)  [psu]
         &       T0   = 40d0, &  ! (in)  [degree C]
         &       Pref = 0d0   &  ! (in)  [dbar]
         &      ) &              ! (out) [degree C]
         &   ) / answer_InSituTemp2Ptemp1 &
         & )

    call  AssertLessThan( &
         &  message="Check SeaWaterProp_FM83_Temp2PTemp", &
         &  answer=acceptableError, check= abs( answer_InSituTemp2Ptemp2 - &
         &     SeaWaterProp_FM83_Temp2PTemp( &
         &       P0   = 1d04, &  ! (in)  [dbar]
         &       S0   = 35d0, &  ! (in)  [psu]
         &       T0   =  0d0, &  ! (in)  [degree C]
         &       Pref =  0d0  &  ! (in)  [dbar]
         &      ) &              ! (out) [degree C]
         &   ) / answer_InSituTemp2Ptemp2 &
         & )
    call  AssertLessThan( &
         &  message="Check SeaWaterProp_FM83_HeatCapacity", &
         &  answer=acceptableError, check= abs( ansert_HeatCapacity - &
         &     SeaWaterProp_FM83_HeatCapacity( &
         &       P   = 1d04, &   ! (in)  [dbar]
         &       S   = 40d0, &   ! (in)  [psu]
         &       T   = 40d0 &    ! (in)  [degree C]
         &      ) &              ! (out) [degree C]
         &   ) / ansert_HeatCapacity &
         & )

    call SeaWaterProp_FM83_Final()

  end subroutine test_FM83


end program SeaWaterProp_mod_test
