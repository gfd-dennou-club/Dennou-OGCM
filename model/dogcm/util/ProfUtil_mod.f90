!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief A module to provide the information of computational time. 
!! 
!! @author Yuta Kawai
!! 
!!
module ProfUtil_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !  
  implicit none
  private
  
  ! 公開手続き
  ! Public procedure
  !
  public :: ProfUtil_Init, ProfUtil_Final
  public :: ProfUtil_RapStart, ProfUtil_RapEnd
  public :: ProfUtil_RapReport

  ! 公開変数
  ! Public variables
  !
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variables
  !
  
  type, public :: ProfUtilInfo_
     character(TOKEN) :: rapName
     real(DP)         :: rapTstr
     real(DP)         :: rapTtot
     integer          :: rapNstr
     integer          :: rapNend
     integer          :: level
  end type ProfUtilInfo_
  
  integer, parameter                :: rapNLimit         = 100
  type(ProfUtilInfo_), save, target :: profUtilInfo(rapNLimit)
  integer                           :: rapNMax           = 0
  integer, parameter                :: DEFAULT_RAP_LEVEL = 2
  integer                           :: rap_level
  
  character(*), parameter:: module_name = 'ProfUtil_mod'
  
contains

  !>
  !!
  !!
  subroutine ProfUtil_Init(configNmlName)

    ! 引用文; Use statement
    !

    !* gtool5
    use dc_iounit, only: FileOpen

    use dc_types, only: STDOUT

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml  
    integer:: iostat_nml
    namelist /ProfUtil_nml/ rap_level
    
    ! 実行文; Executable statements

    ! Default values settings
    !
    rap_level = DEFAULT_RAP_LEVEL

    ! Input from NAMELIST
    !
    if ( trim(configNmlName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlName, mode = 'r' )  ! (in)

       rewind( unit_nml )
       read( unit_nml, &               ! (in)
            & nml = profutil_nml, &    ! (out)
            & iostat = iostat_nml )    ! (out)
       close( unit_nml )
   end if

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  rap_level          = %d', i = (/ rap_level /) )
    
  end subroutine ProfUtil_Init

  !-------------------------------------------------------------------------------------

  !>
  !!
  !!
  subroutine ProfUtil_Final()

    ! 実行文; Executable statements
    !

    rapNMax = 0
    
  end subroutine ProfUtil_Final

  !-------------------------------------------------------------------------------------

  !>
  !!
  !!
  subroutine ProfUtil_RapStart( rapname_base, level )

    ! 引用文; Use statement
    !
    use omp_lib
    
    ! 宣言文; Declaration statement
    !
    character(*), intent(in)      :: rapname_base
    integer, intent(in), optional :: level
    
    ! 局所変数
    ! Local variables
    !
    integer :: id
    integer :: level_
    character(TOKEN) :: rapname
    type(ProfUtilInfo_), pointer :: prfInfo
    
    ! 実行文; Executable statements
    !

    if (present(level)) then
       level_ = level
    else
       level_ = DEFAULT_RAP_LEVEL
    end if
    
    if ( level_ > rap_level ) return

    rapname = rapname_base
    id = get_rapid( rapname, level_ )

    prfInfo => profUtilInfo(id)    
    prfInfo%rapTstr = omp_get_wtime()
    prfInfo%rapNstr = prfInfo%rapNstr + 1

  end subroutine ProfUtil_RapStart
  
  !>
  !!
  !!
  subroutine ProfUtil_RapEnd( rapname_base, level )

    ! 引用文; Use statement
    !
    use omp_lib
    
    ! 宣言文; Declaration statement
    !
    character(*), intent(in)      :: rapname_base
    integer, intent(in), optional :: level
    
    ! 局所変数
    ! Local variables
    !
    integer :: id
    integer :: level_
    character(TOKEN) :: rapname
    type(ProfUtilInfo_), pointer :: prfInfo
    
    ! 実行文; Executable statements
    !

    if (present(level)) then
       if ( level > rap_level ) return
    end if
    
    rapname = rapname_base
    id = get_rapid( rapname, level_ )
    if (level_ > rap_level) return

    prfInfo => profUtilInfo(id)
    prfInfo%rapTtot = prfInfo%rapTtot + (omp_get_wtime() -  prfInfo%rapTstr)
    prfInfo%rapNend = prfInfo%rapNend + 1
    
  end subroutine ProfUtil_RapEnd

  !>
  !!
  !!
  subroutine ProfUtil_RapReport()

    ! 引用文; Use statement
    !
    
    !* gtool5
    use dc_string, only: &
         & CPrintf
    
    use dc_iounit, only: FileOpen
    
    ! 局所変数
    ! Local variables
    !
    integer :: id
    type(ProfUtilInfo_), pointer :: prfInfo

    character(*), parameter :: PROF_FNAME = "./Dennou-OGCM.prof"
    character(STRING) :: buff
    integer :: unit00
    
    ! 実行文; Executable statements
    !

    do id=1, rapNMax
       prfInfo => profUtilInfo(id)
       if ( prfInfo%rapNstr /= prfInfo%rapNend ) then
          call MessageNotify('E', module_name, "Mismatch report: rapname=%a", &
               & ca=(/ prfInfo%rapName /) )
       end if
    end do

    ! Output the information of computational time
    !
    call MessageNotify('M', module_name, "Output the information of computational time to '%a'.", &
         & ca=(/ PROF_FNAME  /) )
    
    call FileOpen( unit00, file=PROF_FNAME, mode='w' )

    write(unit00,*) "= ProfUtil setting"
    buff = CPrintf("rap_level=%d, time_units='sec'", i=(/ rap_level /))
    write(unit00,*) trim(buff)
    write(unit00,*)

    write(unit00,*) "= Computational time"    
    do id=1, rapNMax
       prfInfo => profUtilInfo(id)       
       buff = CPrintf( "[%a] Time=%r (N=%d)", &
            & ca=(/ prfInfo%rapName /), r=(/ real(prfInfo%rapTtot) /), i=(/ prfInfo%rapNstr /) )
       write(unit00,*) trim(buff)
    end do

    close( unit00 )
    
  end subroutine ProfUtil_RapReport
  
  !- Private subroutines ---------------------------------------------------------------

  function get_rapid( rapname, level ) result(id)

    ! 宣言文; Declaration statement
    !    
    character(*), intent(in)    :: rapname
    integer,      intent(inout) :: level
    integer                     :: id

    ! 局所変数
    ! Local variables
    !
    character(TOKEN) :: trapname
    
    ! 実行文; Executable statements
    !

    trapname = trim(rapname)
    
    do id=1, rapNMax
       if( trapname == profUtilInfo(id)%rapname ) then
          level = profUtilInfo(id)%level
          return
       end if
    end do

    rapNMax = rapNMax + 1
    id      = rapNMax
    profUtilInfo(id)%rapName = trapname
    profUtilInfo(id)%rapNstr = 0
    profUtilInfo(id)%rapNend = 0
    profUtilInfo(id)%rapTtot = 0d0
    profUtilInfo(id)%level   = level

  end function get_rapid
  
end module ProfUtil_mod

