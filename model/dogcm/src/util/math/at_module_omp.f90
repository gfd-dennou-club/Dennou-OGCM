!--
!----------------------------------------------------------------------
! Copyright(c) 2002-2010 SPMDODEL Development Group. All rights reserved.
!----------------------------------------------------------------------
!
!ɽ��  at_module
!
!      spml/at_module �⥸�塼��� 1 ����ͭ���ΰ�β��Ǥ�ή�α�ư��
!      �����ӥ������Ѵ��ˤ�ꥹ�ڥ��ȥ���ͷ׻����뤿��� Fortran90 �ؿ���
!      �󶡤���.
!
!      2 �����ǡ����� 1 �����˴ؤ���Ʊ���˥��ڥ��ȥ�׻���¹Ԥ��뤿���
!      �ؿ����󶡤��Ƥ���, 2, 3 �����ΰ�Ǥη׻��Υ١������󶡤���.
!
!      ���Υ⥸�塼��������� ISPACK/FTPACK �� Fortran77 ���֥롼�����
!      �Ƥ�Ǥ���. ���ڥ��ȥ�ǡ�������ӳʻ����ǡ����γ�Ǽ��ˡ�ˤĤ��Ƥ�
!      ISPACK/FTPACK �Υޥ˥奢��򻲾Ȥ��줿��.
!
!
!����  2002/01/24  �ݹ�����
!      2002/02/06  �ݹ�����  �Ť�Ƴ��. ���̥��󥿡��ե������ʤ���.
!                            ����η��˱��������󥿡��ե�����̾
!      2002/03/25  �ݹ�����  �⥸�塼��̾�ѹ�
!      2002/08/20  �ݹ�����  ��ʬ��ʿ�Ѵؿ��ɲ�
!      2002/11/10  �ݹ�����  �¶��֤Ǥζ������롼�����ɲ�
!      2005/01/09  �ݹ�����  msgdmp -> MessageNotify ���ѹ�
!      2006/03/06  �ݹ�����  �����Ȥ� RDoc �Ѥ˽���
!      2006/03/19  �ݹ�����  �ѿ�����³����������򥳥��Ȥ��ɲ�
!      2007/10/24  �ݹ�����  ����Ѵؿ��ɲ�
!                            at_Initial �Ƥ�ľ�����Ȳ�ǽ�ˤ��뤿���
!                            allocate ���줿 public �ѿ��� deallocate ����.
!      2007/11/21  �ݹ�����  ��������֥롼�����å���������
!      2007/12/27  �ݹ�����  �����ȴְ㤤����
!      2009/01/04  �ݹ�����  spml ��ˡ��ȿ���뤿��
!                            Interpolate1dim_t, a_Interpolate1dim_at ���
!      2009/01/09  �ݹ�����  at_Initial ��å����������դ��ɲ�
!      2009/01/29  ��������ʿ �����Ȥ�RDoc �Ѥ˽���
!      2009/07/31  �ݹ�����  �������׻�������� threadprivate ����(OpenMP)
!      2009/12/05  �ݹ�����  threadprivate ���ꥳ���ȥ�����
!      2010/03/10  ��������ʿ  threadprivate ���(����ѥ����¸)
!
!++
module at_module_omp
  !
  != at_module
  !
  ! Authors:: Shin-ichi Takehiro, Youhei SASAKI
  ! Version:: $Id: at_module.f90 598 2013-08-20 03:23:44Z takepiro $
  ! Copyright&License:: See COPYRIGHT[link:../COPYRIGHT]
  !
  !== ����
  !
  ! spml/at_module �⥸�塼��� 1 ����ͭ���ΰ�β��Ǥ�ή�α�ư��
  ! �����ӥ������Ѵ��ˤ�ꥹ�ڥ��ȥ���ͷ׻����뤿��� Fortran90 �ؿ���
  ! �󶡤���.
  !
  ! 2 �����ǡ����� 1 �����˴ؤ���Ʊ���˥��ڥ��ȥ�׻���¹Ԥ��뤿���
  ! �ؿ����󶡤��Ƥ���, 2, 3 �����ΰ�Ǥη׻��Υ١������󶡤���.
  !
  ! ���Υ⥸�塼��������� ISPACK/FTPACK �� Fortran77 ���֥롼�����
  ! �Ƥ�Ǥ���. ���ڥ��ȥ�ǡ�������ӳʻ����ǡ����γ�Ǽ��ˡ�ˤĤ��Ƥ�
  ! ISPACK/FTPACK �Υޥ˥奢��򻲾Ȥ��줿��.
  !
  !
  !== �ؿ����ѿ���̾���ȷ��ˤĤ���
  !
  !=== ̿̾ˡ
  !
  ! * �ؿ�̾����Ƭ (t_, g_, at_, ag_) ��, �֤��ͤη��򼨤��Ƥ���.
  !   t_  :: �����ӥ����եǡ���
  !   g_  :: 1 �����ʻ����ǡ���
  !   at_ :: 1 ���������ӥ����եǡ�����ʣ���¤�� 2 �����ǡ���
  !   ag_ :: 1 �����ʻ����ǡ�����ʣ���¤�� 2 �����ǡ���.
  !
  ! * �ؿ�̾�δ֤�ʸ����(Dx)��, ���δؿ��κ��Ѥ�ɽ���Ƥ���.
  !
  ! * �ؿ�̾�κǸ� (_e,_at,_g, _ag) ��, �����ѿ��η��������ӥ����եǡ���
  !   ����ӳʻ����ǡ����Ǥ��뤳�Ȥ򼨤��Ƥ���.
  !   _t  :: �����ӥ����եǡ���
  !   _g  :: 1 �����ʻ����ǡ���
  !   _at :: 1 ���������ӥ����եǡ�����ʣ���¤�� 2 �����ǡ���
  !   _ag :: 1 �����ʻ����ǡ�����ʣ���¤�� 2 �����ǡ���
  !
  !=== �ƥǡ����μ��������
  !
  ! * g : 1 �����ʻ����ǡ���.
  !   * �ѿ��μ���ȼ����� real(8), dimension(0:im).
  !   * im �� X ��ɸ�γʻ������Ǥ���, ���֥롼���� at_Initial �ˤ�
  !     ���餫�������ꤷ�Ƥ���.
  !
  ! * t : �����ӥ����եǡ���.
  !   * �ѿ��μ���ȼ����� real(8), dimension(0:km).
  !   * km �� X �����κ����ȿ��Ǥ���, ���֥롼���� at_Initial �ˤ�
  !     ���餫�������ꤷ�Ƥ���. ���ڥ��ȥ�ǡ����γ�Ǽ�Τ������ˤĤ��Ƥ�...
  !
  ! * ag : 1 ����(X)�ʻ����ǡ������¤�� 2 �����ǡ���.
  !   * �ѿ��μ���ȼ����� real(8), dimension(:,0:im).
  !     �� 2 ������ X ������ɽ��.
  !
  ! * at : 1 ���������ӥ����եǡ������¤�� 2 �����ǡ���.
  !   * �ѿ��μ���ȼ����� real(8), dimension(:,0:km).
  !     �� 2 ���������ڥ��ȥ��ɽ��.
  !
  ! * g_ �ǻϤޤ�ؿ����֤��ͤ� 1 �����ʻ����ǡ�����Ʊ��.
  !
  ! * t_ �ǻϤޤ�ؿ����֤��ͤϥ����ӥ����եǡ�����Ʊ��.
  !
  ! * ag_ �ǻϤޤ�ؿ����֤��ͤ� 1 �����ʻ����ǡ������¤��
  !   2 �����ǡ�����Ʊ��.
  !
  ! * at_ �ǻϤޤ�ؿ����֤��ͤ� 1 ���������ӥ����եǡ������¤��
  !   2 �����ǡ�����Ʊ��.
  !
  ! * �����ӥ����եǡ������Ф�����ʬ���κ��ѤȤ�, �б�����ʻ����ǡ�����
  !   ��ʬ�ʤɤ���Ѥ������ǡ���������ӥ������Ѵ�������ΤΤ��ȤǤ���.
  !
  !== �ѿ�����³����������
  !
  !==== �����
  !
  ! at_Initial  :: �����ӥ������Ѵ��γʻ�����, �ȿ�, �ΰ���礭��������
  !
  !==== ��ɸ�ѿ�
  !
  ! g_X        :: �ʻ�����ɸ(X)���Ǽ���� 1 ��������
  ! g_X_Weight :: �Ťߺ�ɸ���Ǽ���� 1 ��������
  !
  !==== �����Ѵ�
  !
  ! g_t, ag_at :: �����ӥ����եǡ�������ʻҥǡ����ؤ��Ѵ�
  ! t_g, at_ag :: �ʻҥǡ�����������ӥ����եǡ����ؤ��Ѵ�
  !
  !==== ��ʬ
  !
  ! t_Dx_t, at_Dx_at :: �����ӥ����եǡ����� X ��ʬ����Ѥ�����
  !
  !==== ���
  !
  ! Interpolate_t,
  ! a_Interpolate_at  ::
  ! �����ӥ����եǡ�������Ǥ�դ����Ǥ��ͤ����
  !
  !==== ����������
  !
  ! at_Boundaries_DD, at_Boundaries_DN,
  ! at_Boundaries_ND, at_Boundaries_NN  ::
  ! �ǥ��ꥯ��,�Υ��ޥ󶭳�����Ŭ��
  !
  ! at_BoundariesTau_DD, at_BoundariesTau_DN,
  ! at_BoundariesTau_ND, at_BoundariesTau_NN  ::
  ! �ǥ��ꥯ��, �Υ��ޥ󶭳�����Ŭ��(����ˡ)
  !
  ! at_BoundariesGrid_DD, at_BoundariesGrid_DN,
  ! at_BoundariesGrid_ND, at_BoundariesGrid_NN  ::
  ! �ǥ��ꥯ��, �Υ��ޥ󶭳�����Ŭ��(����ˡ)
  !
  !==== ��ʬ��ʿ��
  !
  ! a_Int_ag, a_Avr_ag :: 1 �����ʻ����ǡ������¤�� 2 �����������ʬ�����ʿ��
  ! Int_g, Avr_g       :: 1 �����ʻ����ǡ�������ʬ�����ʿ��
  !
  use dc_message
  use lumatrix

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none
  private
  public g_X, g_X_Weight                      ! ��ɸ�ѿ�
  public at_Initial                           ! �����
  public ag_at, at_ag, g_t, t_g               ! �����Ѵ�
  public at_Dx_at, t_Dx_t                     ! ��ʬ
  public Interpolate_t, a_Interpolate_at      ! ���
  public at_Boundaries_DD, at_Boundaries_DN   ! �������
  public at_Boundaries_ND, at_Boundaries_NN   ! �������
  public at_BoundariesTau_DD, at_BoundariesTau_DN     ! �������
  public at_BoundariesTau_ND, at_BoundariesTau_NN     ! �������
  public at_BoundariesGrid_DD, at_BoundariesGrid_DN   ! �������
  public at_BoundariesGrid_ND, at_BoundariesGrid_NN   ! �������
  public a_Int_ag, Int_g, a_Avr_ag, Avr_g             ! ��ʬ��ʿ��

  interface at_Boundaries_DD
     !
     ! ξü�ǥ��ꥯ�췿��������Ŭ��(����ˡ).
     ! ξ�����Ǥ��ͤ�Ϳ����.
     !
     ! * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !   �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !   �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     !== �����ȷ�̤η�
     !
     ! * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     ! * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesTau_DD_1d, at_BoundariesTau_DD_2d
  end interface

  interface at_Boundaries_DN
     !
     ! �ǥ��ꥯ�졦�Υ��ޥ󷿶�������Ŭ��(����ˡ).
     ! i=0 ����, i=im �Ǹ��ۤ��ͤ�Ϳ����.
     !
     ! * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !   �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !   �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     !== �����ȷ�̤η�
     !
     ! * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     ! * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesTau_DN_1d, at_BoundariesTau_DN_2d
  end interface

  interface at_Boundaries_ND
     !
     ! �Υ��ޥ󡦥ǥ��ꥯ�췿��������Ŭ��(����ˡ).
     ! i=0 �Ǹ��ۤ���, i=im ���ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     !== �����ȷ�̤η�
     !
     ! * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !   real(8), dimension(:,0:km),intent(inout)       :: at_data
     !   !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !   real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !   !(in) Ŭ�Ѥ��붭����
     !
     ! * 1 ���������ӥ����եǡ����ξ��
     !
     !   real(8), dimension(0:km),intent(inout)       :: t_data
     !   !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !   real(8), dimension(2), intent(in), optional  :: values
     !   !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesTau_ND_1d, at_BoundariesTau_ND_2d
  end interface

  interface at_Boundaries_NN
     !
     ! ξü�Υ��ޥ����Ŭ��(����ˡ).
     ! ξ�����Ǹ��ۤ��ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     !==  �����ȷ�̤η�
     !
     ! * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !   real(8), dimension(:,0:km),intent(inout)       :: at_data
     !   !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !   real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !   !(in) Ŭ�Ѥ��붭����
     !
     ! * 1 ���������ӥ����եǡ����ξ��
     !
     !   real(8), dimension(0:km),intent(inout)       :: t_data
     !   !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !   real(8), dimension(2), intent(in), optional  :: values
     !   !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesTau_NN_1d, at_BoundariesTau_NN_2d
  end interface

  interface at_BoundariesTau_DD
     !
     ! ξü�ǥ��ꥯ�췿��������Ŭ��(����ˡ).
     ! ξ�����Ǥ��ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     !== �����ȷ�̤η�
     !
     !  * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     !  * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesTau_DD_1d, at_BoundariesTau_DD_2d
  end interface

  interface at_BoundariesTau_DN
     !
     ! �ǥ��ꥯ�졦�Υ��ޥ󷿶�������Ŭ��(����ˡ).
     ! i=0 ����, i=im �Ǹ��ۤ��ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     ! �����ȷ�̤η�
     !
     !  * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     !  * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesTau_DN_1d, at_BoundariesTau_DN_2d
  end interface

  interface at_BoundariesTau_ND
     !
     ! �Υ��ޥ󡦥ǥ��ꥯ�췿��������Ŭ��(����ˡ).
     ! i=0 �Ǹ��ۤ���, i=im ���ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     ! �����ȷ�̤η�
     !
     !  * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     !  * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesTau_ND_1d, at_BoundariesTau_ND_2d
  end interface

  interface at_BoundariesTau_NN
     !
     ! ξü�Υ��ޥ����Ŭ��(����ˡ).
     ! ξ�����Ǹ��ۤ��ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     ! �����ȷ�̤η�
     !
     !  * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     !  * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesTau_NN_1d, at_BoundariesTau_NN_2d
  end interface

  interface at_BoundariesGrid_DD
     !
     ! ξü�ǥ��ꥯ�췿��������Ŭ��(�¶��֤Ǥ�ɾ��).
     ! ξ�����Ǥ��ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     ! �����ȷ�̤η�
     !
     !  * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     !  * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesGrid_DD_1d, at_BoundariesGrid_DD_2d
  end interface

  interface at_BoundariesGrid_DN
     !
     ! �ǥ��ꥯ�졦�Υ��ޥ󷿶�������Ŭ��(�¶��֤Ǥ�ɾ��).
     ! i=0 ����, i=im �Ǹ��ۤ��ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     ! �����ȷ�̤η�
     !
     !  * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     !  * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesGrid_DN_1d, at_BoundariesGrid_DN_2d
  end interface

  interface at_BoundariesGrid_ND
     !
     ! �Υ��ޥ󡦥ǥ��ꥯ�췿��������Ŭ��(�¶��֤Ǥ�ɾ��)
     ! i=0 �Ǹ��ۤ���, i=im ���ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     ! �����ȷ�̤η�
     !
     !  * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     !  * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesGrid_ND_1d, at_BoundariesGrid_ND_2d
  end interface

  interface at_BoundariesGrid_NN
     !
     ! ξü�Υ��ޥ����Ŭ��(�¶��֤Ǥ�ɾ��).
     ! ξ�����Ǹ��ۤ��ͤ�Ϳ����.
     !
     !  * ��������Ŭ�Ѥ�������μ����ˤ�ä������ǥ��֥롼�����
     !    �Ȥ�ʬ���Ƥ���. �桼�������󥿡��ե������϶��̤Ǥ���Τ�
     !    �����롼�����Ƥ�ɬ�פϤʤ�.
     !
     ! �����ȷ�̤η�
     !
     !  * 1 ���������ӥ����եǡ������¤�� 2 ��������ξ��
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     !  * 1 ���������ӥ����եǡ����ξ��
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) Ŭ�Ѥ��붭����
     !
     module procedure at_BoundariesGrid_NN_1d, at_BoundariesGrid_NN_2d
  end interface


  integer, dimension(5)              :: it
  real(8), dimension(:), allocatable :: t
  real(8), parameter                 :: pi=3.1415926535897932385D0

  integer :: im, km                        ! �ʻ�����, �����ȿ�
  real(8) :: xl                            ! �ΰ���礭��

  real(8), allocatable :: g_X(:)
  ! �ʻ�����ɸ
  ! km ���Υ����ӥ�����¿�༰������������ޤ�ʻ���

  real(8), allocatable :: g_X_Weight(:)
  ! �ʻ����Ťߺ�ɸ
  ! �Ƴʻ����ˤ�������ʬ�Τ���νŤߤ���Ǽ���Ƥ���

  save :: im, km, xl, it, t, g_X, g_X_Weight
  
  integer :: nThread
  logical :: omp_flag

contains

! ---- ����� ----
  subroutine at_Initial(i,k,xmin,xmax, threadNum)
    !
    ! �����ӥ������Ѵ��γʻ�����, �ȿ�, �ΰ���礭�������ꤹ��.
    !
    ! ¾�δؿ����ѿ���Ƥ�����, �ǽ�ˤ��Υ��֥롼�����Ƥ��
    ! �������򤷤ʤ���Фʤ�ʤ�.
    !
    integer,intent(in) :: i              !(in) �ʻ�����
    integer,intent(in) :: k              !(in) �����ȿ�
    real(8),intent(in) :: xmin, xmax     !(in) ��ɸ���ϰ�
    integer, optional, intent(in) :: threadNum

    integer :: ii,kk
    
    im=i ; km=k
    xl = xmax-xmin

    if ( im <= 0 .or. km <= 0 ) then
       call MessageNotify('E','at_initial', &
            'Number of grid points and waves should be positive')
    elseif ( mod(im,2) /= 0 ) then
       call MessageNotify('E','at_initial','Number of grid points should be even')
    elseif ( km > im ) then
       call MessageNotify('E','at_initial','KM shoud be less equal IM')
    endif

    if ( allocated(t) ) deallocate(t)
    allocate(t(3*im))
    call fttcti(im,it,t)

    if ( allocated(g_X) ) deallocate(g_X)
    allocate(g_X(0:im))
    do ii=0,im
       g_X(ii) = (xmax+xmin)/2 + xl/2 * cos(pi*ii/im)
    enddo

    if ( allocated(g_X_Weight) ) deallocate(g_X_Weight)
    allocate(g_X_Weight(0:im))
    do ii=0,im
       g_X_Weight(ii) = 1.0
       do kk=2,km,2
          g_X_Weight(ii) = g_X_Weight(ii) &
                          + 2/(1D0-kk**2) * cos(kk*ii*pi/im)
       enddo
       if ( (km == im) .and. (mod(im,2)==0) ) then  ! �Ǹ���¤� factor 1/2.
          g_X_Weight(ii) = g_X_Weight(ii) &
                          - 1/(1D0-km**2)* cos(km*ii*pi/im)
       endif
       g_X_Weight(ii) = 2D0/im * g_X_Weight(ii) * xl/2
    enddo
    g_X_Weight(0)  = g_X_Weight(0) / 2
    g_X_Weight(im) = g_X_Weight(im) / 2


    !
    
#ifdef _OPENMP

    omp_flag = .true.
    !$omp parallel
    !$omp single
    nThread = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    if(present(threadNum)) then
       nThread = threadNum
    end if
    call MessageNotify('M', 'at_initial', "Number of thread is set %d.", i=(/nThread/))
#else
    nThread = 1
    omp_flag = .false.
#endif


    call MessageNotify('M','at_initial','at_module (2009/07/31) is initialized')

  end subroutine at_Initial

! ---- ���Ѵ� ----
  function ag_at(at_data)

    use omp_lib
    !
    ! �����ӥ����եǡ�������ʻҥǡ������Ѵ�����(2 ����������).
    !
    real(8), dimension(:,:), intent(in)      :: at_data
    !(in) �����ӥ����եǡ���

    real(8), dimension(size(at_data,1),0:im) :: ag_at
    !(out) �ʻ����ǡ���

    real(8), dimension(size(at_data,1)*im) :: y
    ! ���������
    integer :: m

    integer :: tr
    integer, allocatable :: lc_m(:), lbnd(:), ubnd(:)
    integer :: lc_mm
    real(8), allocatable :: lc_ag_at(:,:), lc_y(:)

    m = size(at_data,1)
    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','ag_at', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','ag_at', &
            'The Chebyshev dimension of input data too large.')
    endif


    if(omp_flag .and. (m > nThread)) then
       allocate( lbnd(nThread), ubnd(nThread), lc_m(nThread) )
       do tr=1, nThread
         lbnd(tr) = (tr-1)*ceiling(m/real(nThread)) + 1
         ubnd(tr) =  min(tr*ceiling(m/real(nThread)), m)
         lc_m(tr) = ubnd(tr) - lbnd(tr) + 1
       end do
       lc_mm = maxval(lc_m)
       allocate( lc_ag_at(lc_mm,0:im), lc_y(lc_mm*im) )
       !$omp parallel do private(lc_ag_at, lc_y)
       do tr=1, nThread
         lc_ag_at(1:lc_m(tr),0:km) = at_data(lbnd(tr):ubnd(tr),:)
         call fttctb(lc_mm,im,lc_ag_at(:,:),lc_y(:),it,t)  
         ag_at(lbnd(tr):ubnd(tr),0:km) = lc_ag_at(lbnd(tr):ubnd(tr),0:km)
       end do
      !  !$omp parallel private(lc_m, lbnd, ubnd, tr)
      !  tr = omp_get_thread_num() + 1
      !  lbnd = (tr-1)*ceiling(m/real(nThread)) + 1
      !  ubnd = min(tr*ceiling(m/real(nThread)), m)
      !  lc_m = ubnd - lbnd + 1

      !  ag_at(lbnd:ubnd,:) = 0d0
      !  ag_at(lbnd:ubnd,0:km) = at_data(lbnd:ubnd,:)
      !  call fttctb(lc_m,im,ag_at(lbnd:ubnd,:),y(im*(lbnd-1)+1:im*ubnd),it,t)
      !  !$omp end parallel
    else
       ag_at = 0.0D0
       ag_at(:,0:km) = at_data
       call fttctb(m,im,ag_at,y,it,t)
    end if
    

  end function ag_at

  function g_t(t_data)
    !
    ! �����ӥ����եǡ�������ʻҥǡ������Ѵ�����(1 ����������).
    !
    real(8), dimension(:), intent(in)  :: t_data
    !(in) �����ӥ����եǡ���

    real(8), dimension(0:im)           :: g_t
    !(out) �ʻ����ǡ���

    real(8), dimension(1,size(t_data)) :: t_work
    ! ���������
    real(8), dimension(1,0:im)         :: g_work
    ! ���������

    t_work(1,:) = t_data
    g_work = ag_at(t_work)
    g_t = g_work(1,:)

  end function g_t


! ---- ���Ѵ� ----
  function at_ag(ag_data)
    !
    ! �ʻҥǡ�����������ӥ����եǡ������Ѵ�����(2 ����������).
    !
    real(8), dimension(:,:), intent(in)      :: ag_data
    !(in) �ʻ����ǡ���

    real(8), dimension(size(ag_data,1),0:km) :: at_ag
    !(out) �����ӥ����եǡ���

    real(8), dimension(size(ag_data,1)*im)   :: y
    real(8), dimension(size(ag_data,1),0:im) :: ag_work
    integer :: m

    integer :: tr
    integer, allocatable :: lc_m(:), lbnd(:), ubnd(:)
    integer :: lc_mm
    real(8), allocatable :: lc_ag_work(:,:), lc_y(:)


    m = size(ag_data,1)
    if ( size(ag_data,2)-1 < im ) then
       call MessageNotify('E','at_ag', &
            'The Grid points of input data too small.')
    elseif ( size(ag_data,2)-1 > im ) then
       call MessageNotify('W','at_ag', &
            'The Grid points of input data too large.')
    endif

    if( omp_flag .and. (m > nThread) ) then
       allocate( lbnd(nThread), ubnd(nThread), lc_m(nThread) )
       do tr=1, nThread
         lbnd(tr) = (tr-1)*ceiling(m/real(nThread)) + 1
         ubnd(tr) =  min(tr*ceiling(m/real(nThread)), m)
         lc_m(tr) = ubnd(tr) - lbnd(tr) + 1
       end do
       lc_mm = maxval(lc_m)
       allocate( lc_ag_work(lc_mm,0:im), lc_y(lc_mm*im) )
       !$omp parallel do private(lc_ag_work, lc_y)
       do tr=1, nThread
         lc_ag_work(1:lc_m(tr),:) = ag_data(lbnd(tr):ubnd(tr),:)
         call fttctf(lc_mm,im,lc_ag_work(:,:),lc_y(:),it,t)  
         at_ag(lbnd(tr):ubnd(tr),:) = lc_ag_work(lbnd(tr):ubnd(tr),0:km)
       end do      
      !  !$omp parallel private(lc_m, lbnd, ubnd, tr)
      !  tr = omp_get_thread_num() + 1
      !  lbnd = (tr-1)*ceiling(m/real(nThread)) + 1
      !  ubnd = min(tr*ceiling(m/real(nThread)), m)
      !  lc_m = ubnd - lbnd + 1

      !  ag_work(lbnd:ubnd,:) = ag_data(lbnd:ubnd,:) 
      !  call fttctf(lc_m,im,ag_work(lbnd:ubnd,:),y(im*(lbnd-1)+1:im*ubnd),it,t)
      !  at_ag(lbnd:ubnd,:) = ag_work(lbnd:ubnd, 0:km)
      !  !$omp end parallel
    else
       ag_work = ag_data
       call fttctf(m,im,ag_work,y,it,t)
       at_ag = ag_work(:,0:km)
    end if

  end function at_ag

  function t_g(g_data)  ! ����ʻ� -> ���ڥ��ȥ�
    !
    ! �ʻҥǡ�����������ӥ����եǡ������Ѵ�����(1 ����������).
    !
    real(8), dimension(:), intent(in)     :: g_data
    !(in) �ʻ����ǡ���

    real(8), dimension(0:km)              :: t_g
    !(out) �����ӥ����եǡ���

    real(8), dimension(1,size(g_data)) :: ag_work
    real(8), dimension(1,0:km)         :: at_work

    ag_work(1,:) = g_data
    at_work = at_ag(ag_work)
    t_g = at_work(1,:)

  end function t_g

! ---- ��ʬ�׻� ----

  function at_Dx_at(at_data_) result(at_Dx)
    !
    ! ���ϥ����ӥ����եǡ����� X ��ʬ����Ѥ���(2 ����������).
    !
    ! �����ӥ����եǡ����� X ��ʬ�Ȥ�, �б�����ʻ����ǡ����� X ��ʬ��
    ! ���Ѥ������ǡ����Υ����ӥ������Ѵ��Τ��ȤǤ���.
    !
    !
    real(8), dimension(:,0:), intent(in)                    :: at_data_
    !(in) ���ϥ����ӥ����եǡ���

    real(8), dimension(size(at_data_,1),0:size(at_data_,2)-1) :: at_Dx
    !(out) �����ӥ����եǡ����� X ��ʬ

    integer :: m, k
    integer :: nm, kmax

    integer :: tr, lc_m, lbnd, ubnd

    nm=size(at_data_,1)
    kmax=size(at_data_,2)-1
    if ( kmax  < km ) then
       call MessageNotify('W','at_Dx_at', &
            'The Chebyshev dimension of input data too small.')
    elseif ( kmax > km ) then
       call MessageNotify('E','at_Dx_at', &
            'The Chebyshev dimension of input data too large.')
    endif


    if(omp_flag .and. (size(at_data_,1) > nThread)) then
!!$       !$omp parallel private(lbnd, ubnd, tr)
!!$       tr = omp_get_thread_num() + 1
!!$       lbnd = (tr-1)*ceiling(nm/real(nThread)) + 1
!!$       ubnd = min(tr*ceiling(nm/real(nThread)), nm)
!!$       at_Dx(lbnd:ubnd,:) = at_Dx_at_core(at_data_(lbnd:ubnd,:))
!!$       !$omp end parallel

       !$omp parallel do 
       do m=1,nm
          at_Dx(m:m,:) = at_Dx_at_core(at_data_(m:m,:))
       end do
    else
       at_Dx(:,:) = at_Dx_at_core(at_data_)
    end if
    

  contains
    function at_Dx_at_core(at_data) result(at_Dx_at)
      !
      ! ���ϥ����ӥ����եǡ����� X ��ʬ����Ѥ���(2 ����������).
      !
      ! �����ӥ����եǡ����� X ��ʬ�Ȥ�, �б�����ʻ����ǡ����� X ��ʬ��
      ! ���Ѥ������ǡ����Υ����ӥ������Ѵ��Τ��ȤǤ���.
      !
      !
      real(8), dimension(:,0:), intent(in)                    :: at_data
      !(in) ���ϥ����ӥ����եǡ���
      
      real(8), dimension(size(at_data,1),0:size(at_data,2)-1) :: at_Dx_at
      !(out) �����ӥ����եǡ����� X ��ʬ

      integer :: m, k
      integer :: nm, kmax

      nm=size(at_data,1)
      kmax=size(at_data,2)-1

      if ( kmax == im ) then
         at_Dx_at(:,kmax)   = 0.
         at_Dx_at(:,kmax-1) = 2 * km * at_data(:,kmax) /2
      else
         at_Dx_at(:,kmax)   = 0.
         at_Dx_at(:,kmax-1) = 2 * km * at_data(:,kmax)
         ! �������Ȥϥ���å��б������ȿ�̤��. Factor 1/2 ����
      endif

      do k=kmax-2,0,-1
         at_Dx_at(:,k) = at_Dx_at(:,k+2) + 2*(k+1)*at_data(:,k+1)
      enddo

      at_Dx_at = 2d0/xl * at_Dx_at

  end function at_Dx_at_core
  end function at_Dx_at

  function t_Dx_t(t_data)
    !
    ! ���ϥ����ӥ����եǡ����� X ��ʬ����Ѥ���(1 ����������).
    !
    ! �����ӥ����եǡ����� X ��ʬ�Ȥ�, �б�����ʻ����ǡ����� X ��ʬ��
    ! ���Ѥ������ǡ����Υ����ӥ������Ѵ��Τ��ȤǤ���.
    !
    !
    real(8), dimension(:), intent(in)   :: t_data
    !(in) ���ϥ����ӥ����եǡ���

    real(8), dimension(size(t_data))    :: t_Dx_t
    !(out) �����ӥ����եǡ����� X ��ʬ

    real(8), dimension(1,size(t_data))  :: at_work
    ! ���������

    at_work(1,:) = t_data
    at_work = at_Dx_at(at_work)
    t_Dx_t = at_work(1,:)

  end function t_Dx_t

! ---- ��ַ׻� ----
  function Interpolate_t(t_data,xval)
    real(8), dimension(0:), intent(in)  :: t_data
    !(in) ���ϥ����ӥ����եǡ���

    real(8), intent(in)                 :: xval
    ! ��֤������κ�ɸ

    real(8)                             :: Interpolate_t
    ! ��֤�����̤���

    integer :: kmax
    ! ��������κ��缡��

    real(8) :: y2, y1, y0, x
    ! Crenshow's reccurence formula �׻����ѿ�

    integer :: k
    ! DO ʸ�ѿ�

    kmax = size(t_data)-1

    x =(xval -(g_X(0)+g_X(im))/2 )/(g_X(0)-g_X(im))*2

    y2 = 0 ; y1 = 0
    do k=kmax,1,-1
       y0 = 2*x*y1 - y2 + t_data(k)
       y2 = y1 ; y1 = y0
    enddo

    Interpolate_t = - y2 + x*y1 + t_data(0)/2
    if ( kmax == im ) then
       Interpolate_t = Interpolate_t -t_data(kmax)/2*cos(kmax*acos(x))
    endif

  end function Interpolate_t

  function a_Interpolate_at(at_data,xval)
    real(8), dimension(:,0:), intent(in) :: at_data
    !(in) ���ϥ����ӥ����եǡ���

    real(8), intent(in)                  :: xval
    ! ��֤������κ�ɸ

    real(8), dimension(size(at_data,1))  :: a_Interpolate_at
    ! ��֤�����̤���

    integer :: kmax
    ! ��������κ��缡��

    real(8), dimension(size(at_data,1))  :: y2, y1, y0
    real(8)                              ::  x
    ! Crenshow's reccurence formula �׻����ѿ�

    integer :: k
    ! DO ʸ�ѿ�

    kmax = size(at_data,2)-1

    x =(xval -(g_X(0)+g_X(im))/2 )/(g_X(0)-g_X(im))*2

    y2 = 0 ; y1 = 0
    do k=kmax,1,-1
       y0 = 2*x*y1 - y2 + at_data(:,k)
       y2 = y1 ; y1 = y0
    enddo

    a_Interpolate_at = - y2 + x*y1 + at_data(:,0)/2
    if ( kmax == im ) then
       a_Interpolate_at = a_Interpolate_at -at_data(:,kmax)/2*cos(kmax*acos(x))
    endif
  end function a_Interpolate_at

!---- Dirichlet ���������(����ˡ) ----

  subroutine at_BoundariesTau_DD_2d(at_data,values)
    !
    ! Dirichlet ����������Ŭ��(����ˡ, 2 ����������)
    ! ξ�����Ǥ��ͤ�Ϳ����.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data ! �ǡ���(m,0:km)
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) ������(m,2)

    real(8), dimension(:,:), allocatable  :: alu
    integer, dimension(:), allocatable    :: kp
    real(8), dimension(0:km,0:km)         :: tt_data
    real(8), dimension(0:km,0:im)         :: tg_data
    real(8), dimension(size(at_data,1))    :: value1, value2           ! ������

    logical :: first = .true.
    integer :: k
    save    :: alu, kp, first

    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','at_BoundariesTau_DD', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','at_BoundariesTau_DD', &
            'The Chebyshev dimension of input data too large.')
    endif

    if (.not. present(values)) then
       value1=0 ; value2=0
    else
       value1 = values(:,1) ; value2 = values(:,2)
    endif

    if ( first ) then
       first = .false.

       allocate(alu(0:km,0:km),kp(0:km))

       tt_data=0
       do k=0,km
          tt_data(k,k)=1
       enddo
       alu = tt_data

       tg_data = ag_at(tt_data)
       alu(km-1,:) = tg_data(:,0)
       alu(km,:)   = tg_data(:,im)

       call ludecomp(alu,kp)
    endif

    at_data(:,km-1) = value1
    at_data(:,km)   = value2
    at_data = lusolve(alu,kp,at_data)

  end subroutine at_BoundariesTau_DD_2d

  subroutine at_BoundariesTau_DD_1d(t_data,values)
    !
    ! Dirichlet ����������Ŭ��(����ˡ, 1 ����������)
    ! ξ�����Ǥ��ͤ�Ϳ����.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) ������

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! ������

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesTau_DD_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesTau_DD_1d


!---- Dirichlet/Neumann ���������(����ˡ) ----

  subroutine at_BoundariesTau_DN_2d(at_data,values)
    !
    ! Dirichlet/Neumann ����������Ŭ��(����ˡ, 2 ����������)
    ! i=0 ����, i=im �Ǹ��ۤ��ͤ�Ϳ����.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) ������(m,2)

    real(8), dimension(:,:), allocatable  :: alu
    integer, dimension(:), allocatable    :: kp
    real(8), dimension(0:km,0:km)         :: tt_data
    real(8), dimension(0:km,0:im)         :: tg_data
    real(8), dimension(size(at_data,1))    :: value1, value2           ! ������

    logical :: first = .true.
    integer :: k
    save    :: alu, kp, first

    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','at_BoundariesTau_DN', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','at_BoundariesTau_DN', &
            'The Chebyshev dimension of input data too large.')
    endif

    if (.not. present(values)) then
       value1=0 ; value2=0
    else
       value1 = values(:,1) ; value2 = values(:,2)
    endif

    if ( first ) then
       first = .false.

       allocate(alu(0:km,0:km),kp(0:km))

       tt_data = 0
       do k=0,km
          tt_data(k,k)=1
       enddo
       alu = tt_data

       tg_data = ag_at(tt_data)
       alu(km-1,:) = tg_data(:,0)
       tg_data = ag_at(at_Dx_at(tt_data))
       alu(km,:)   = tg_data(:,im)

       call ludecomp(alu,kp)
    endif

    at_data(:,km-1) = value1
    at_data(:,km)   = value2
    at_data = lusolve(alu,kp,at_data)

  end subroutine at_BoundariesTau_DN_2d

  subroutine at_BoundariesTau_DN_1d(t_data,values)
    !
    ! Dirichlet/Neumann ����������Ŭ��(����ˡ, 1 ����������)
    ! i=0 ����, i=im �Ǹ��ۤ��ͤ�Ϳ����.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) ������

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! ������

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesTau_DN_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesTau_DN_1d

!---- Neumann/Dirichlet ���������(����ˡ) ----

  subroutine at_BoundariesTau_ND_2d(at_data,values)
    !
    ! Neumann/Dirichlet ����������Ŭ��(����ˡ, 2 ����������)
    ! i=0 �Ǹ��ۤ���, i=im ���ͤ�Ϳ����.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) ������(m,2)

    real(8), dimension(:,:), allocatable  :: alu
    integer, dimension(:), allocatable    :: kp
    real(8), dimension(0:km,0:km)         :: tt_data
    real(8), dimension(0:km,0:im)         :: tg_data
    real(8), dimension(size(at_data,1))    :: value1, value2           ! ������

    logical :: first = .true.
    integer :: k
    save    :: alu, kp, first

    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','at_BoundariesTau_ND', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','at_BoundariesTau_ND', &
            'The Chebyshev dimension of input data too large.')
    endif

    if (.not. present(values)) then
       value1=0 ; value2=0
    else
       value1 = values(:,1) ; value2 = values(:,2)
    endif

    if ( first ) then
       first = .false.

       allocate(alu(0:km,0:km),kp(0:km))

       tt_data = 0
       do k=0,km
          tt_data(k,k)=1
       enddo
       alu = tt_data

       tg_data = ag_at(at_Dx_at(tt_data))
       alu(km-1,:) = tg_data(:,0)
       tg_data = ag_at(tt_data)
       alu(km,:)   = tg_data(:,im)

       call ludecomp(alu,kp)
    endif

    at_data(:,km-1) = value1
    at_data(:,km)   = value2
    at_data = lusolve(alu,kp,at_data)

  end subroutine at_BoundariesTau_ND_2d

  subroutine at_BoundariesTau_ND_1d(t_data,values)
    !
    ! Neumann/Dirichlet ����������Ŭ��(����ˡ, 1 ����������)
    ! i=0 �Ǹ��ۤ���, i=im ���ͤ�Ϳ����.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) ������

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! ������

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesTau_ND_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesTau_ND_1d


!---- Neumann ���������(����ˡ) ----

  subroutine at_BoundariesTau_NN_2d(at_data,values)
    !
    ! Neumann/Dirichlet ����������Ŭ��(����ˡ, 2 ����������)
    ! ξ�����Ǹ��ۤ��ͤ�Ϳ����.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) ������(m,2)

    real(8), dimension(:,:), allocatable  :: alu
    integer, dimension(:), allocatable    :: kp
    real(8), dimension(0:km,0:km)         :: tt_data
    real(8), dimension(0:km,0:im)         :: tg_data
    real(8), dimension(size(at_data,1))    :: value1, value2           ! ������

    logical :: first = .true.
    integer :: k
    save    :: alu, kp, first

    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','at_BoundariesTau_NN', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','at_BoundariesTau_NN', &
            'The Chebyshev dimension of input data too large.')
    endif

    if (.not. present(values)) then
       value1=0 ; value2=0
    else
       value1 = values(:,1) ; value2 = values(:,2)
    endif

    if ( first ) then
       first = .false.

       allocate(alu(0:km,0:km),kp(0:km))

       tt_data = 0
       do k=0,km
          tt_data(k,k)=1
       enddo
       alu = tt_data

       tg_data = ag_at(at_Dx_at(tt_data))
       alu(km-1,:) = tg_data(:,0)
       alu(km,:)   = tg_data(:,im)

       call ludecomp(alu,kp)
    endif

    at_data(:,km-1) = value1
    at_data(:,km)   = value2
    at_data = lusolve(alu,kp,at_data)

  end subroutine at_BoundariesTau_NN_2d

  subroutine at_BoundariesTau_NN_1d(t_data,values)
    !
    ! Neumann/Dirichlet ����������Ŭ��(����ˡ, 1 ����������)
    ! ξ�����Ǹ��ۤ��ͤ�Ϳ����.
    !
    ! ���Υ��֥롼�����ľ�ܻȤ����Ȥ򴫤�ʤ�.
    ! ���̥��󥿡��ե����� at_Boundaries_NN ���Ѥ��뤳��.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) ������

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! ������

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesTau_NN_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesTau_NN_1d

  !--------------- ��ʬ�׻� -----------------
    function a_Int_ag(ag)
      !
      ! 1 �����ʻ����ǡ������¤�� 2 �����������ʬ
      !
      real(8), dimension(:,0:), intent(in)     :: ag
      !(in)���ϳʻ����ǡ���

      real(8), dimension(size(ag,1))           :: a_Int_ag
      !(out) ��ʬ�����ǡ���
      integer :: i

      if ( size(ag,2) < im+1 ) then
         call MessageNotify('E','ae_ag', &
              'The Grid points of input data too small.')
      elseif ( size(ag,2) > im+1 ) then
         call MessageNotify('W','ae_ag', &
              'The Grid points of input data too large.')
      endif

      a_Int_ag = 0.0d0
      do i=0,im
         a_Int_ag(:) = a_Int_ag(:) + ag(:,i)*g_X_Weight(i)
      enddo
    end function a_Int_ag

    function Int_g(g)
      !
      ! 1 �����ʻ����ǡ�������ʬ�����ʿ��.
      !
      real(8), dimension(0:im), intent(in)   :: g
      !(in) �ʻ����ǡ���

      real(8)                                :: Int_g
      !(out) ��ʬ��

      Int_g = sum(g*g_X_Weight)
    end function Int_g

    function a_Avr_ag(ag)
      !
      ! 1 �����ʻ����ǡ������¤�� 2 ���������ʿ��
      !
      real(8), dimension(:,0:), intent(in)   :: ag
      !(in)���ϳʻ����ǡ���

      real(8), dimension(size(ag,1))         :: a_Avr_ag
      !(out) ʿ�Ѥ����ǡ���

      a_Avr_ag = a_Int_ag(ag)/sum(g_X_Weight)
    end function a_Avr_ag

    function Avr_g(g)
      !
      ! 1 �����ʻ����ǡ�����ʿ��
      !
      real(8), dimension(0:im), intent(in)   :: g
      !(in) �ʻ����ǡ���

      real(8)                                :: Avr_g
      !(out) ��ʬ��

      Avr_g = Int_g(g)/sum(g_X_Weight)
    end function Avr_g



!---- Dirichlet ���������(�¶��֤Ǥ�ɾ��) ----

  subroutine at_BoundariesGrid_DD_2d(at_data,values)
    !
    ! Dirichlet ����������Ŭ��(�¶��֤Ǥ�ɾ��, 2 ����������)
    ! ξ�����Ǥ��ͤ�Ϳ����.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) ������(m,2)

    real(8), dimension(:,:), allocatable     :: alu
    integer, dimension(:), allocatable       :: kp
    real(8), dimension(size(at_data,1),0:im) :: ag_data
    real(8), dimension(0:km,0:km)            :: tt_data
    real(8), dimension(0:km,0:im)            :: tg_data
    real(8), dimension(size(at_data,1))      :: value1, value2  ! ������

    logical :: first = .true.
    integer :: k
    save    :: alu, kp, first

    if ( im /= km ) then
       call MessageNotify('E','at_BoundariesGrid_DD', &
            'Chebyshev truncation and number of grid points should be same.')
    endif

    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','at_BoundariesGrid_DD', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','at_BoundariesGrid_DD', &
            'The Chebyshev dimension of input data too large.')
    endif

    if (.not. present(values)) then
       value1=0 ; value2=0
    else
       value1 = values(:,1) ; value2 = values(:,2)
    endif

    if ( first ) then
       first = .false.
       allocate(alu(0:im,0:km),kp(0:im))

       tt_data = 0
       do k=0,km
          tt_data(k,k)=1.0
       enddo
       tg_data = ag_at(tt_data)
       alu = transpose(tg_data)
!       alu(km-1,:) = tg_data(:,0)
!       alu(km,:)   = tg_data(:,im)

       call ludecomp(alu,kp)
    endif

    ag_data = ag_at(at_data)
    ag_data(:,0)  = value1
    ag_data(:,im) = value2
    at_data = lusolve(alu,kp,ag_data)

  end subroutine at_BoundariesGrid_DD_2d

  subroutine at_BoundariesGrid_DD_1d(t_data,values)
    !
    ! Dirichlet ����������Ŭ��(�¶��֤Ǥ�ɾ��, 1 ����������)
    ! ξ�����Ǥ��ͤ�Ϳ����.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) ������

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! ������

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesGrid_DD_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesGrid_DD_1d

!---- Dirichlet/Neumann ���������(�¶��֤Ǥ�ɾ��) ----

  subroutine at_BoundariesGrid_DN_2d(at_data,values)
    !
    ! Dirichlet/Neumann ����������Ŭ��(�¶��֤Ǥ�ɾ��, 2 ����������)
    ! i=0 ����, i=im �Ǹ��ۤ��ͤ�Ϳ����.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) ������(m,2)

    real(8), dimension(:,:), allocatable     :: alu
    integer, dimension(:), allocatable       :: kp
    real(8), dimension(size(at_data,1),0:im) :: ag_data
    real(8), dimension(0:km,0:km)            :: tt_data
    real(8), dimension(0:km,0:im)            :: tg_data
    real(8), dimension(size(at_data,1))      :: value1, value2  ! ������

    logical :: first = .true.
    integer :: k
    save    :: alu, kp, first

    if ( im /= km ) then
       call MessageNotify('E','at_BoundariesGrid_DN', &
            'Chebyshev truncation and number of grid points should be same.')
    endif

    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','at_BoundariesGrid_DN', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','at_BoundariesGrid_DN', &
            'The Chebyshev dimension of input data too large.')
    endif

    if (.not. present(values)) then
       value1=0 ; value2=0
    else
       value1 = values(:,1) ; value2 = values(:,2)
    endif

    if ( first ) then
       first = .false.
       allocate(alu(0:im,0:km),kp(0:im))

       tt_data = 0
       do k=0,km
          tt_data(k,k)=1.0
       enddo
       tg_data = ag_at(tt_data)
       alu = transpose(tg_data)

       tg_data = ag_at(at_dx_at(tt_data))
       alu(im,:) = tg_data(:,im)

       call ludecomp(alu,kp)
    endif

    ag_data = ag_at(at_data)
    ag_data(:,0)  = value1
    ag_data(:,im) = value2
    at_data = lusolve(alu,kp,ag_data)

  end subroutine at_BoundariesGrid_DN_2d

  subroutine at_BoundariesGrid_DN_1d(t_data,values)
    !
    ! Dirichlet/Neumann ����������Ŭ��(�¶��֤Ǥ�ɾ��, 1 ����������)
    ! i=0 ����, i=im �Ǹ��ۤ��ͤ�Ϳ����.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) ������

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! ������

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesGrid_DN_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesGrid_DN_1d

!---- Neumann/Dirichlet ���������(�¶��֤Ǥ�ɾ��) ----

  subroutine at_BoundariesGrid_ND_2d(at_data,values)
    !
    ! Neumann/Dirichlet ����������Ŭ��(�¶��֤Ǥ�ɾ��, 2 ����������)
    ! i=0 �Ǹ��ۤ���, i=im ���ͤ�Ϳ����.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) ������(m,2)

    real(8), dimension(:,:), allocatable     :: alu
    integer, dimension(:), allocatable       :: kp
    real(8), dimension(size(at_data,1),0:im) :: ag_data
    real(8), dimension(0:km,0:km)            :: tt_data
    real(8), dimension(0:km,0:im)            :: tg_data
    real(8), dimension(size(at_data,1))      :: value1, value2  ! ������

    logical :: first = .true.
    integer :: k
    save    :: alu, kp, first

    if ( im /= km ) then
       call MessageNotify('E','at_BoundariesGrid_ND', &
            'Chebyshev truncation and number of grid points should be same.')
    endif

    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','at_BoundariesGrid_ND', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','at_BoundariesGrid_DD', &
            'The Chebyshev dimension of input data too large.')
    endif

    if (.not. present(values)) then
       value1=0 ; value2=0
    else
       value1 = values(:,1) ; value2 = values(:,2)
    endif

    if ( first ) then
       first = .false.
       allocate(alu(0:im,0:km),kp(0:im))

       tt_data = 0
       do k=0,km
          tt_data(k,k)=1.0
       enddo
       tg_data = ag_at(tt_data)
       alu = transpose(tg_data)

       tg_data = ag_at(at_dx_at(tt_data))
       alu(0,:)  = tg_data(:,0)

       call ludecomp(alu,kp)
    endif

    ag_data = ag_at(at_data)
    ag_data(:,0)  = value1
    ag_data(:,im) = value2
    at_data = lusolve(alu,kp,ag_data)

  end subroutine at_BoundariesGrid_ND_2d

  subroutine at_BoundariesGrid_ND_1d(t_data,values)
    !
    ! Neumann ����������Ŭ��(�¶��֤Ǥ�ɾ��, 1 ����������)
    ! i=0 �Ǹ��ۤ���, i=im ���ͤ�Ϳ����.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) ������

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! ������

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesGrid_ND_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesGrid_ND_1d

!---- Neumann ��������� ----

  subroutine at_BoundariesGrid_NN_2d(at_data,values)
    !
    ! Neumann ����������Ŭ��(�¶��֤Ǥ�ɾ��, 2 ����������)
    ! ξ�����Ǹ��ۤ��ͤ�Ϳ����.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) ������(m,2)

    real(8), dimension(:,:), allocatable     :: alu
    integer, dimension(:), allocatable       :: kp
    real(8), dimension(size(at_data,1),0:im) :: ag_data
    real(8), dimension(0:km,0:km)            :: tt_data
    real(8), dimension(0:km,0:im)            :: tg_data
    real(8), dimension(size(at_data,1))      :: value1, value2  ! ������

    logical :: first = .true.
    integer :: k
    save    :: alu, kp, first

    if ( im /= km ) then
       call MessageNotify('E','at_BoundariesGrid_NN', &
            'Chebyshev truncation and number of grid points should be same.')
    endif

    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','at_BoundariesGrid_NN', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','at_BoundariesGrid_NN', &
            'The Chebyshev dimension of input data too large.')
    endif

    if (.not. present(values)) then
       value1=0 ; value2=0
    else
       value1 = values(:,1) ; value2 = values(:,2)
    endif

    if ( first ) then
       first = .false.
       allocate(alu(0:im,0:km),kp(0:im))

       tt_data = 0
       do k=0,km
          tt_data(k,k)=1.0
       enddo
       tg_data = ag_at(tt_data)
       alu = transpose(tg_data)

       tg_data = ag_at(at_dx_at(tt_data))
       alu(0,:)  = tg_data(:,0)
       alu(im,:) = tg_data(:,im)

       call ludecomp(alu,kp)
    endif

    ag_data = ag_at(at_data)
    ag_data(:,0)  = value1
    ag_data(:,im) = value2
    at_data = lusolve(alu,kp,ag_data)

  end subroutine at_BoundariesGrid_NN_2d

  subroutine at_BoundariesGrid_NN_1d(t_data,values)
    !
    ! Neumann ����������Ŭ��(�¶��֤Ǥ�ɾ��, 1 ����������)
    ! ξ�����Ǹ��ۤ��ͤ�Ϳ����.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) ��������Ŭ�Ѥ�������ӥ����եǡ���(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) ������

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! ������

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesGrid_NN_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesGrid_NN_1d

end module at_module_omp
