module OutputData_mod

  use dc_types
  use dc_message
  use gtool_history

  use VectorSpace_mod
  use SphericalCoord_mod
  use PolyMesh_mod
  use GeometricField_mod
  use fvMeshInfo_mod
  use fvCalculus_mod
  use netcdfDataWriter_mod

  use SimParameters_mod

  use GridSet_mod

  use VariableSet_mod, only: &
       & v_h, v_hb, s_normalVel, v_div

  implicit none

  public :: OutputData, OutputDataAnalysisInfo

  type(netcdfDataWriter), save :: ncWriter
  type(PointScalarField), save :: p_zeta

contains

subroutine OutputData_Init()

  !
  !
  call GeometricField_Init(p_zeta, plMesh, "p_zeta", "relative vorticity", "s-1")

  call netcdfDataWriter_Init(ncWriter, "data.nc", plMesh)
  call netcdfDataWriter_Regist( ncWriter, (/ v_h, v_div /), (/ p_zeta /) )

  !
  !
  !
  call HistoryCreate( &               
    & file='dataAnalysis.nc', title='data analysis', &
    & source='Sample program of gtool_history/gtool5',   &
    & institution='GFD_Dennou Club davis project',       &
    & dims=(/'time'/), dimsizes=(/0/),               &
    & longnames=(/'time'/),       &
    & units=(/'s'/),                                 &
    & origin=real(0), interval=real(outputIntrVal) )


  call HistoryAddVariable( &
    & varname='l2norm', dims=(/'time'/), &
    & longname='l2norm for analystic solution', units='1', xtype='double')

  call HistoryAddVariable( &
    & varname='linfnorm', dims=(/'time'/), &
    & longname='linfnorm for analystic solution', units='1', xtype='double')

  call HistoryAddVariable( &
    & varname='KEDblTime', dims=(/'time'/), &
    & longname='KE doubling time scale', units='year', xtype='double')


end subroutine OutputData_Init


subroutine OutputData_Final()

  call netcdfDataWriter_Final(ncWriter)

  call GeometricField_Final(p_zeta)

end subroutine OutputData_Final

subroutine OutputData( tstep )

  integer, intent(in) :: tstep
  character(STRING) :: dataFileName

integer :: maxId(2)


maxId = maxloc(v_h%data%v_)
write(*,*) "max height=", At(v_h,maxId(2)), ":", RadToDegUnit(CartToSphPos(plmesh%cellPosList(maxId(2))))


  p_zeta = curl(s_normalVel)
  v_div = div(s_normalVel)
  
  call netcdfDataWriter_write(ncWriter, v_h)
  call netcdfDataWriter_write(ncWriter, v_div)
  call netcdfDataWriter_write(ncWriter, p_zeta)
  call netcdfDataWriter_AdvanceTimeStep(ncWriter, dble(tstep*delTime) )

end subroutine OutputData

subroutine OutputDataAnalysisInfo(tstep, & 
  & l2norm, linfnorm )

  integer, intent(in) :: tstep
  real(DP), optional :: l2norm, linfnorm

  if( present(l2norm) ) call HistoryPut("l2norm", l2norm)
  if( present(linfnorm) ) call HistoryPut("linfnorm", linfnorm)
 
end subroutine OutputDataAnalysisInfo

end module OutputData_mod
