program SpmlUtil_mod_test
  
  use dc_types
  use dc_message
  use Constants_mod
  use GridSet_mod
  use SpmlUtil_mod
  use dc_test
  use dc_string

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

#ifdef DSOGCM_MODE_AXISYM
  character(*), parameter :: configNmlFile = "defaultConfig_axisym.nml"
#else
  character(*), parameter :: configNmlFile = "defaultConfig.nml"
#endif


  integer :: caseId

  ! Variables for continious field
  integer :: nMode
  real(DP), parameter :: contField_bt_ErrLims(8) = (/ 1d-15, 1d-15, 1d-15, 1d-13, 1d-13, 1d-9, 1d-8, 1d-7 /) 
  real(DP), parameter :: contField_st_ErrLims(8) = (/ 1d-14, 1d-14, 1d-14, 1d-12, 1d-10, 1d-9, 1d-8, 1d-7 /) 

  ! Variables for noncontinious field  
  real(DP), parameter :: lyrLenRatios(5) = (/ 0.05d0, 0.01d0, 0.005d0, 0.0025d0, 0.001d0 /)
  real(DP), parameter :: ncontField_bt_ErrLims(5) = (/ 1d-11, 1d-5, 1d-3, 1d-3, 1d-3 /)
  real(DP), parameter :: ncontField_st_ErrLims(5) = (/ 1d-12, 1d-5, 1d-4, 1d-3, 1d-3 /)

  real(DP) :: lyrLen
  integer :: nThread

  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  call Constants_Init(configNmlFile)
  call GridSet_Init(configNmlFile)

#ifdef _OPENMP
  !$omp parallel
  !$omp single
  nThread = omp_get_num_threads()
  !$omp end single
  !$omp end parallel 
#else
  nThread = 1
#endif

  call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet, nThread)
  call GridSet_construct()

  call MessageNotify("M", "SpmlUtil_mod_test", "= WaveFunc")
  do caseId=1, size(contField_bt_ErrLims)
     nMode = caseId
     call check_xy_IntSig_BtmToTop_xyz( &
          & eval_func1_intBtmToTop, check_intfunc1_intBtmToTop, caseId, contField_bt_ErrLims(caseId) )
  end do
  do caseId=1, 8
     nMode = caseId
     call check_xyz_IntSig_SigToTop_xyz( &
          & eval_func1_intSigToTop, check_intfunc1_intSigToTop, caseId, CPrintf("waveFunc_%d_k%d.dat",i=(/nMode,kMax/)), &
          & contField_st_ErrLims(caseId) )
  end do

  call MessageNotify("M", "SpmlUtil_mod_test", "= ExpFunc")
  do caseId=1, size(lyrLenRatios)
     lyrLen = lyrLenRatios(caseId)
     call check_xy_IntSig_BtmToTop_xyz( &
          & eval_func2_intBtmToTop, check_intfunc2_intBtmToTop, caseId, ncontField_bt_ErrLims(caseId) )
  end do
  do caseId=1,size(lyrLenRatios)

     lyrLen = lyrLenRatios(caseId)
     call check_xyz_IntSig_SigToTop_xyz( &
          & eval_func2_intSigToTop, check_intfunc2_intSigToTop, caseId, CPrintf("expFunc_%d_k%d.dat",i=(/caseId,kMax/)), &
          & ncontField_st_ErrLims(caseId) )
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
    val = sin(nMode*PI)/(nMode*PI)
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
    val = exp(sig/lyrLen) + exp(-(1d0+sig)/lyrLen)
  end function eval_func2_intBtmToTop
  
  pure function check_intfunc2_intBtmToTop() result(val)
    real(DP) :: val
    val = 2d0*(1d0-exp(-1d0/lyrLen))*lyrLen
  end function check_intfunc2_intBtmToTop

  pure function eval_func2_intSigToTop(sig) result(val)
    real(DP), intent(in) :: sig
    real(DP) :: val

    val = exp(sig/lyrLen) + exp(-(1d0+sig)/lyrLen)
  end function eval_func2_intSigToTop
  
  pure function check_intfunc2_intSigToTop(sig) result(val)
    real(DP), intent(in) :: sig
    real(DP) :: val

    val = lyrLen*(1d0 - exp(sig/lyrLen) - exp(-1d0/lyrLen) + exp(-(1d0+sig)/lyrLen) )
  end function check_intfunc2_intSigToTop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

  subroutine check_xy_IntSig_BtmToTop_xyz(eval_func, check_intfunc, caseId, errorLimit)

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

    integer, intent(in) :: caseId
    real(DP), intent(in) :: errorLimit

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
    message=CPrintf("xy_IntSig_BtmToTop_xyz: param=%d, l2norm=%f", i=(/caseId/), d=(/ l2norm /))
    call  AssertLessThan(message=message, &
         & answer=errorLimit, check=l2norm)

  end subroutine check_xy_IntSig_BtmToTop_xyz

  subroutine check_xyz_IntSig_SigToTop_xyz(eval_func, check_intfunc, caseId, outputFileName, errorLimit)

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
    integer, intent(in) :: caseId
    character(*), optional :: outputFileName
    real(DP), intent(in) :: errorLimit

    real(DP) :: val(0:iMax-1,jMax,0:kMax)
    real(DP) :: intVal(0:iMax-1,jMax, 0:kMax)
    real(DP) :: checkVal(0:iMax-1,jMax, 0:kMax)
    real(DP) :: lostCheck(0:iMax-1,jMax, 0:kMax)
    integer :: k, m
    real(DP) :: l2norm, totIntVal
    character(STRING) :: message

    do k=0, kMax
       val(:,:,k) = eval_func(g_Sig(k))
       checkVal(:,:,k) = check_intfunc(g_Sig(k))
    end do

    intVal = xyz_IntSig_SigToTop_xyz(val)
    l2norm = abs(sum(intVal - checkVal))/dble(iMax*jMax*(kMax+1))       
    message=CPrintf("xyz_IntSig_SigToTop_xyz: param %d, l2norm %f", i=(/caseId/), d=(/ l2norm/))
    lostCheck = (val - xyz_wt(-wt_DSig_wt(wt_xyz(intVal))))/maxval(abs(val))
    
    call  AssertLessThan(message=message, answer=errorLimit, check=l2norm)

    !
    
    if (present(outputFileName)) then
       call MessageNotify("M", "SpmlUtil_mod_test", "Output error data..")
       open(10, file=trim(outputFileName), status="replace")
       
       write(10,*) "# Sig,   f(Sig),  normalized integralation value,  nomalized analytic integration value, l2norm of error,", &
            & "information loss with spectral conversion", "Zero"

       do k=0,kMax
          write(10,'(7E15.5e4)') g_Sig(k), eval_func(g_Sig(k)), &
               & intVal(1,1,k)/maxval(abs(checkVal(1,1,:))), checkVal(1,1,k)/maxval(abs(checkVal(1,1,:))), & 
               & sqrt( (intVal(1,1,k) - checkVal(1,1,k))**2 )/maxval(abs(checkVal(1,1,:))), lostCheck(1,1,k), 0d0
       end do
       close(10)
    end if



  end subroutine check_xyz_IntSig_SigToTop_xyz

end program SpmlUtil_mod_test
