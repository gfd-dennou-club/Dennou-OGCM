program SpmlUtil_mod_test
  
  use dc_types
  use dc_message
  use Constants_mod
  use GridSet_mod
  use SpmlUtil_mod
  use dc_test
  use dc_string

  implicit none

  character(*), parameter :: configNmlFile = "defaultConfig.nml"
  integer :: nMode
  integer :: lyrDivId
  real(DP) :: lyrLen
  real(DP) :: lyrLenRatio(6) = (/ 2d0, 1d0, 0.5d0, 0.25d0, 0.1d0, 0.05d0 /)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  call Constants_Init(configNmlFile)
  call GridSet_Init(configNmlFile)
  call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
  call GridSet_construct()

  call MessageNotify("M", "SpmlUtil_mod_test", "= WaveFunc")
  do nMode=1, 5
     call check_xy_IntSig_BtmToTop_xyz( &
          & eval_func1_intBtmToTop, check_intfunc1_intBtmToTop, dble(nMode) )
  end do
  do nMode=1, 8
     call check_xyz_IntSig_SigToTop_xyz( &
          & eval_func1_intSigToTop, check_intfunc1_intSigToTop, dble(nMode), CPrintf("waveFunc_%d.dat",i=(/nMode/)) )
  end do

  call MessageNotify("M", "SpmlUtil_mod_test", "= ExpFunc")
  do lyrDivId=1, size(lyrLenRatio)
     lyrLen = lyrLenRatio(lyrDivId)/dble(kMax)
     call check_xy_IntSig_BtmToTop_xyz( &
          & eval_func2_intBtmToTop, check_intfunc2_intBtmToTop, dble(lyrDivId) )
  end do
  do lyrDivId=1,size(lyrLenRatio)

     lyrLen = lyrLenRatio(lyrDivId)/dble(kMax)
     call check_xyz_IntSig_SigToTop_xyz( &
          & eval_func2_intSigToTop, check_intfunc2_intSigToTop, dble(lyrDivId), CPrintf("expFunc_%d.dat",i=(/lyrDivId/)) )
  end do

  call SpmlUtil_Final()
  call GridSet_Final() 
  call Constants_Final()

contains

  !***************************************
  !* type1 
  !***************************************

  pure function eval_func1_intBtmToTop(sig) result(val)

    real(DP), intent(in) :: sig
    real(DP) :: val

    val = cos(nMode*PI*(1d0 + sig))    
  end function eval_func1_intBtmToTop
  
  pure function check_intfunc1_intBtmToTop() result(val)
    real(DP) :: val
    val = sin(nMode*PI)
  end function check_intfunc1_intBtmToTop

  pure function eval_func1_intSigToTop(sig) result(val)
    real(DP), intent(in) :: sig
    real(DP) :: val

    val = cos(nMode*PI*(1d0 + sig))    
  end function eval_func1_intSigToTop
  
  pure function check_intfunc1_intSigToTop(sig) result(val)
    real(DP), intent(in) :: sig
    real(DP) :: val
    val = -sin(nMode*PI*(1d0 + sig))/(nMode*PI)
  end function check_intfunc1_intSigToTop

  !***************************************
  !* type2 
  !***************************************

  pure function eval_func2_intBtmToTop(sig) result(val)

    real(DP), intent(in) :: sig
    real(DP) :: val
    val = exp(sig/lyrLen) - exp(-(1d0+sig)/lyrLen)
  end function eval_func2_intBtmToTop
  
  pure function check_intfunc2_intBtmToTop() result(val)
    real(DP) :: val
    val = 0d0
  end function check_intfunc2_intBtmToTop

  pure function eval_func2_intSigToTop(sig) result(val)
    real(DP), intent(in) :: sig
    real(DP) :: val

    val = exp(sig/lyrLen) - exp(-(1d0+sig)/lyrLen)
  end function eval_func2_intSigToTop
  
  pure function check_intfunc2_intSigToTop(sig) result(val)
    real(DP), intent(in) :: sig
    real(DP) :: val

    val = lyrLen*(1d0 - exp(sig/lyrLen) + exp(-1d0/lyrLen) - exp(-(1d0+sig)/lyrLen) )
  end function check_intfunc2_intSigToTop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

  subroutine check_xy_IntSig_BtmToTop_xyz(eval_func, check_intfunc, param)

    interface
       pure function eval_func(sig) result(val)
         use dc_types
         real(DP), intent(in) :: sig
         real(DP) :: val
       end function eval_func
       pure function check_intfunc() result(val)
         use dc_types
         real(DP) :: val
       end function check_intfunc
    end interface

    real(DP), intent(in) :: param

    real(DP) :: val(0:iMax-1,jMax,0:kMax)
    integer :: k
    real(DP) :: l2norm
    character(STRING) :: message
    real(DP) :: intVal(0:iMax-1,jMax), ansVal
 
    do k=0, kMax
       val(:,:,k) = eval_func(g_Sig(k))
    end do
    intVal = xy_IntSig_BtmToTop_xyz(val)
    l2norm = abs(intVal(1,1) - check_intfunc())
    message=CPrintf("xy_IntSig_BtmToTop_xyz: param=%f, l2norm=%f", d=(/ param, l2norm /))
    call  AssertLessThan(message=message, &
         & answer=5d-10, check=l2norm)

  end subroutine check_xy_IntSig_BtmToTop_xyz

  subroutine check_xyz_IntSig_SigToTop_xyz(eval_func, check_intfunc, param, outputFileName)

    interface
       pure function eval_func(sig) result(val)
         use dc_types
         real(DP), intent(in) :: sig
         real(DP) :: val
       end function eval_func
       pure function check_intfunc(sig) result(val)
         use dc_types
         real(DP), intent(in) :: sig
         real(DP) :: val
       end function check_intfunc
    end interface
    real(DP), intent(in) :: param
    character(*), optional :: outputFileName

    real(DP) :: val(0:iMax-1,jMax,0:kMax)
    real(DP) :: intVal(0:iMax-1,jMax, 0:kMax)
    real(DP) :: checkVal(0:iMax-1,jMax, 0:kMax)
    integer :: k, m
    real(DP) :: theta(0:kMax), l2norm
    character(STRING) :: message

    theta = PI*(1d0 + g_Sig)

    do k=0, kMax
       val(:,:,k) = eval_func(g_Sig(k))
       checkVal(:,:,k) = check_intfunc(g_Sig(k))
    end do
    intVal = xyz_IntSig_SigToTop_xyz(val)

    l2norm = abs(sum(intVal - checkVal))/dble(iMax*jMax*(kMax+1))       
    message=CPrintf("xyz_IntSig_SigToTop_xyz: param %f, l2norm %f", d=(/param, l2norm/))


!!$do k=0,kMax
!!$   write(*,*) "k=",k, ": theta=", theta(k), "val=", val(10,10,k), ": integral ",  intVal(10,10,k), checkVal(10,10,k)
!!$end do

    !
    
    if (present(outputFileName)) then
       call MessageNotify("M", "SpmlUtil_mod_test", "Output error data..")
       open(10, file=trim(outputFileName), status="replace")
       do k=0,kMax
          write(10,'(5E15.5e4)') g_Sig(k), eval_func(g_Sig(k)), intVal(1,1,k), checkVal(1,1,k), & 
               & sqrt( (intVal(1,1,k) - checkVal(1,1,k))**2 )
       end do
       close(10)
    end if

    !
    call  AssertLessThan(message=message, answer=1d-0, check=l2norm)


  end subroutine check_xyz_IntSig_SigToTop_xyz

end program SpmlUtil_mod_test
