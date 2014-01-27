!--
!----------------------------------------------------------------------
! Copyright(c) 2002-2010 SPMDODEL Development Group. All rights reserved.
!----------------------------------------------------------------------
!
!表題  at_module
!
!      spml/at_module モジュールは 1 次元有限領域の下での流体運動を
!      チェビシェフ変換によりスペクトル数値計算するための Fortran90 関数を
!      提供する.
!
!      2 次元データの 1 次元に関して同時にスペクトル計算を実行するための
!      関数も提供しており, 2, 3 次元領域での計算のベースも提供する.
!
!      このモジュールは内部で ISPACK/FTPACK の Fortran77 サブルーチンを
!      呼んでいる. スペクトルデータおよび格子点データの格納方法については
!      ISPACK/FTPACK のマニュアルを参照されたい.
!
!
!履歴  2002/01/24  竹広真一
!      2002/02/06  竹広真一  重み導入. 共通インターフェースなくす.
!                            配列の形に応じたインターフェース名
!      2002/03/25  竹広真一  モジュール名変更
!      2002/08/20  竹広真一  積分・平均関数追加
!      2002/11/10  竹広真一  実空間での境界条件ルーチン追加
!      2005/01/09  竹広真一  msgdmp -> MessageNotify に変更
!      2006/03/06  竹広真一  コメントを RDoc 用に修正
!      2006/03/19  竹広真一  変数・手続き群の要約をコメントに追加
!      2007/10/24  竹広真一  補間用関数追加
!                            at_Initial 呼び直すこと可能にするために
!                            allocate された public 変数を deallocate する.
!      2007/11/21  竹広真一  初期化サブルーチンメッセージ出力
!      2007/12/27  竹広真一  コメント間違い訂正
!      2009/01/04  竹広真一  spml 書法に反するため
!                            Interpolate1dim_t, a_Interpolate1dim_at 削除
!      2009/01/09  竹広真一  at_Initial メッセージに日付を追加
!      2009/01/29  佐々木洋平 コメントをRDoc 用に修正
!      2009/07/31  竹広真一  境界条件計算用配列を threadprivate 指定(OpenMP)
!      2009/12/05  竹広真一  threadprivate 指定コメントアウト
!      2010/03/10  佐々木洋平  threadprivate 削除(コンパイラ依存)
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
  !== 概要
  !
  ! spml/at_module モジュールは 1 次元有限領域の下での流体運動を
  ! チェビシェフ変換によりスペクトル数値計算するための Fortran90 関数を
  ! 提供する.
  !
  ! 2 次元データの 1 次元に関して同時にスペクトル計算を実行するための
  ! 関数も提供しており, 2, 3 次元領域での計算のベースも提供する.
  !
  ! このモジュールは内部で ISPACK/FTPACK の Fortran77 サブルーチンを
  ! 呼んでいる. スペクトルデータおよび格子点データの格納方法については
  ! ISPACK/FTPACK のマニュアルを参照されたい.
  !
  !
  !== 関数・変数の名前と型について
  !
  !=== 命名法
  !
  ! * 関数名の先頭 (t_, g_, at_, ag_) は, 返す値の形を示している.
  !   t_  :: チェビシェフデータ
  !   g_  :: 1 次元格子点データ
  !   at_ :: 1 次元チェビシェフデータが複数並んだ 2 次元データ
  !   ag_ :: 1 次元格子点データが複数並んだ 2 次元データ.
  !
  ! * 関数名の間の文字列(Dx)は, その関数の作用を表している.
  !
  ! * 関数名の最後 (_e,_at,_g, _ag) は, 入力変数の形がチェビシェフデータ
  !   および格子点データであることを示している.
  !   _t  :: チェビシェフデータ
  !   _g  :: 1 次元格子点データ
  !   _at :: 1 次元チェビシェフデータが複数並んだ 2 次元データ
  !   _ag :: 1 次元格子点データが複数並んだ 2 次元データ
  !
  !=== 各データの種類の説明
  !
  ! * g : 1 次元格子点データ.
  !   * 変数の種類と次元は real(8), dimension(0:im).
  !   * im は X 座標の格子点数であり, サブルーチン at_Initial にて
  !     あらかじめ設定しておく.
  !
  ! * t : チェビシェフデータ.
  !   * 変数の種類と次元は real(8), dimension(0:km).
  !   * km は X 方向の最大波数であり, サブルーチン at_Initial にて
  !     あらかじめ設定しておく. スペクトルデータの格納のされ方については...
  !
  ! * ag : 1 次元(X)格子点データの並んだ 2 次元データ.
  !   * 変数の種類と次元は real(8), dimension(:,0:im).
  !     第 2 次元が X 方向を表す.
  !
  ! * at : 1 次元チェビシェフデータの並んだ 2 次元データ.
  !   * 変数の種類と次元は real(8), dimension(:,0:km).
  !     第 2 次元がスペクトルを表す.
  !
  ! * g_ で始まる関数が返す値は 1 次元格子点データに同じ.
  !
  ! * t_ で始まる関数が返す値はチェビシェフデータに同じ.
  !
  ! * ag_ で始まる関数が返す値は 1 次元格子点データの並んだ
  !   2 次元データに同じ.
  !
  ! * at_ で始まる関数が返す値は 1 次元チェビシェフデータの並んだ
  !   2 次元データに同じ.
  !
  ! * チェビシェフデータに対する微分等の作用とは, 対応する格子点データに
  !   微分などを作用させたデータをチェビシェフ変換したもののことである.
  !
  !== 変数・手続き群の要約
  !
  !==== 初期化
  !
  ! at_Initial  :: チェビシェフ変換の格子点数, 波数, 領域の大きさの設定
  !
  !==== 座標変数
  !
  ! g_X        :: 格子点座標(X)を格納した 1 次元配列
  ! g_X_Weight :: 重み座標を格納した 1 次元配列
  !
  !==== 基本変換
  !
  ! g_t, ag_at :: チェビシェフデータから格子データへの変換
  ! t_g, at_ag :: 格子データからチェビシェフデータへの変換
  !
  !==== 微分
  !
  ! t_Dx_t, at_Dx_at :: チェビシェフデータに X 微分を作用させる
  !
  !==== 補間
  !
  ! Interpolate_t,
  ! a_Interpolate_at  ::
  ! チェビシェフデータから任意の点での値を求める
  !
  !==== 境界値問題
  !
  ! at_Boundaries_DD, at_Boundaries_DN,
  ! at_Boundaries_ND, at_Boundaries_NN  ::
  ! ディリクレ,ノイマン境界条件の適用
  !
  ! at_BoundariesTau_DD, at_BoundariesTau_DN,
  ! at_BoundariesTau_ND, at_BoundariesTau_NN  ::
  ! ディリクレ, ノイマン境界条件の適用(タウ法)
  !
  ! at_BoundariesGrid_DD, at_BoundariesGrid_DN,
  ! at_BoundariesGrid_ND, at_BoundariesGrid_NN  ::
  ! ディリクレ, ノイマン境界条件の適用(選点法)
  !
  !==== 積分・平均
  !
  ! a_Int_ag, a_Avr_ag :: 1 次元格子点データの並んだ 2 次元配列の積分および平均
  ! Int_g, Avr_g       :: 1 次元格子点データの積分および平均
  !
  use dc_message
  use lumatrix
  implicit none
  private
  public g_X, g_X_Weight                      ! 座標変数
  public at_Initial                           ! 初期化
  public ag_at, at_ag, g_t, t_g               ! 基本変換
  public at_Dx_at, t_Dx_t                     ! 微分
  public Interpolate_t, a_Interpolate_at      ! 補間
  public at_Boundaries_DD, at_Boundaries_DN   ! 境界条件
  public at_Boundaries_ND, at_Boundaries_NN   ! 境界条件
  public at_BoundariesTau_DD, at_BoundariesTau_DN     ! 境界条件
  public at_BoundariesTau_ND, at_BoundariesTau_NN     ! 境界条件
  public at_BoundariesGrid_DD, at_BoundariesGrid_DN   ! 境界条件
  public at_BoundariesGrid_ND, at_BoundariesGrid_NN   ! 境界条件
  public a_Int_ag, Int_g, a_Avr_ag, Avr_g             ! 積分・平均

  interface at_Boundaries_DD
     !
     ! 両端ディリクレ型境界条件の適用(タウ法).
     ! 両境界での値を与える.
     !
     ! * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !   使い分けている. ユーザーインターフェースは共通であるので
     !   下部ルーチンを呼ぶ必要はない.
     !
     !== 引数と結果の型
     !
     ! * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     ! * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesTau_DD_1d, at_BoundariesTau_DD_2d
  end interface

  interface at_Boundaries_DN
     !
     ! ディリクレ・ノイマン型境界条件の適用(タウ法).
     ! i=0 で値, i=im で勾配の値を与える.
     !
     ! * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !   使い分けている. ユーザーインターフェースは共通であるので
     !   下部ルーチンを呼ぶ必要はない.
     !
     !== 引数と結果の型
     !
     ! * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     ! * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesTau_DN_1d, at_BoundariesTau_DN_2d
  end interface

  interface at_Boundaries_ND
     !
     ! ノイマン・ディリクレ型境界条件の適用(タウ法).
     ! i=0 で勾配の値, i=im で値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     !== 引数と結果の型
     !
     ! * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !   real(8), dimension(:,0:km),intent(inout)       :: at_data
     !   !(inout) 境界条件を適用するチェビシェフデータ
     !
     !   real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !   !(in) 適用する境界値
     !
     ! * 1 次元チェビシェフデータの場合
     !
     !   real(8), dimension(0:km),intent(inout)       :: t_data
     !   !(inout) 境界条件を適用するチェビシェフデータ
     !
     !   real(8), dimension(2), intent(in), optional  :: values
     !   !(in) 適用する境界値
     !
     module procedure at_BoundariesTau_ND_1d, at_BoundariesTau_ND_2d
  end interface

  interface at_Boundaries_NN
     !
     ! 両端ノイマン条件の適用(タウ法).
     ! 両境界で勾配の値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     !==  引数と結果の型
     !
     ! * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !   real(8), dimension(:,0:km),intent(inout)       :: at_data
     !   !(inout) 境界条件を適用するチェビシェフデータ
     !
     !   real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !   !(in) 適用する境界値
     !
     ! * 1 次元チェビシェフデータの場合
     !
     !   real(8), dimension(0:km),intent(inout)       :: t_data
     !   !(inout) 境界条件を適用するチェビシェフデータ
     !
     !   real(8), dimension(2), intent(in), optional  :: values
     !   !(in) 適用する境界値
     !
     module procedure at_BoundariesTau_NN_1d, at_BoundariesTau_NN_2d
  end interface

  interface at_BoundariesTau_DD
     !
     ! 両端ディリクレ型境界条件の適用(タウ法).
     ! 両境界での値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     !== 引数と結果の型
     !
     !  * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     !  * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesTau_DD_1d, at_BoundariesTau_DD_2d
  end interface

  interface at_BoundariesTau_DN
     !
     ! ディリクレ・ノイマン型境界条件の適用(タウ法).
     ! i=0 で値, i=im で勾配の値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     ! 引数と結果の型
     !
     !  * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     !  * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesTau_DN_1d, at_BoundariesTau_DN_2d
  end interface

  interface at_BoundariesTau_ND
     !
     ! ノイマン・ディリクレ型境界条件の適用(タウ法).
     ! i=0 で勾配の値, i=im で値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     ! 引数と結果の型
     !
     !  * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     !  * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesTau_ND_1d, at_BoundariesTau_ND_2d
  end interface

  interface at_BoundariesTau_NN
     !
     ! 両端ノイマン条件の適用(タウ法).
     ! 両境界で勾配の値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     ! 引数と結果の型
     !
     !  * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     !  * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesTau_NN_1d, at_BoundariesTau_NN_2d
  end interface

  interface at_BoundariesGrid_DD
     !
     ! 両端ディリクレ型境界条件の適用(実空間での評価).
     ! 両境界での値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     ! 引数と結果の型
     !
     !  * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     !  * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesGrid_DD_1d, at_BoundariesGrid_DD_2d
  end interface

  interface at_BoundariesGrid_DN
     !
     ! ディリクレ・ノイマン型境界条件の適用(実空間での評価).
     ! i=0 で値, i=im で勾配の値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     ! 引数と結果の型
     !
     !  * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     !  * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesGrid_DN_1d, at_BoundariesGrid_DN_2d
  end interface

  interface at_BoundariesGrid_ND
     !
     ! ノイマン・ディリクレ型境界条件の適用(実空間での評価)
     ! i=0 で勾配の値, i=im で値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     ! 引数と結果の型
     !
     !  * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     !  * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesGrid_ND_1d, at_BoundariesGrid_ND_2d
  end interface

  interface at_BoundariesGrid_NN
     !
     ! 両端ノイマン条件の適用(実空間での評価).
     ! 両境界で勾配の値を与える.
     !
     !  * 境界条件を適用する配列の次元によって内部でサブルーチンを
     !    使い分けている. ユーザーインターフェースは共通であるので
     !    下部ルーチンを呼ぶ必要はない.
     !
     ! 引数と結果の型
     !
     !  * 1 次元チェビシェフデータの並んだ 2 次元配列の場合
     !
     !    real(8), dimension(:,0:km),intent(inout)       :: at_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(size(at_data,1),2), intent(in), optional :: values
     !    !(in) 適用する境界値
     !
     !  * 1 次元チェビシェフデータの場合
     !
     !    real(8), dimension(0:km),intent(inout)       :: t_data
     !    !(inout) 境界条件を適用するチェビシェフデータ
     !
     !    real(8), dimension(2), intent(in), optional  :: values
     !    !(in) 適用する境界値
     !
     module procedure at_BoundariesGrid_NN_1d, at_BoundariesGrid_NN_2d
  end interface


  integer, dimension(5)              :: it
  real(8), dimension(:), allocatable :: t
  real(8), parameter                 :: pi=3.1415926535897932385D0

  integer :: im, km                        ! 格子点数, 切断波数
  real(8) :: xl                            ! 領域の大きさ

  real(8), allocatable :: g_X(:)
  ! 格子点座標
  ! km 次のチェビシェフ多項式の零点から定まる格子点

  real(8), allocatable :: g_X_Weight(:)
  ! 格子点重み座標
  ! 各格子点における積分のための重みが格納してある

  save :: im, km, xl, it, t, g_X, g_X_Weight

contains

! ---- 初期化 ----
  subroutine at_Initial(i,k,xmin,xmax)
    !
    ! チェビシェフ変換の格子点数, 波数, 領域の大きさを設定する.
    !
    ! 他の関数や変数を呼ぶ前に, 最初にこのサブルーチンを呼んで
    ! 初期設定をしなければならない.
    !
    integer,intent(in) :: i              !(in) 格子点数
    integer,intent(in) :: k              !(in) 切断波数
    real(8),intent(in) :: xmin, xmax     !(in) 座標の範囲

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
       if ( (km == im) .and. (mod(im,2)==0) ) then  ! 最後の和は factor 1/2.
          g_X_Weight(ii) = g_X_Weight(ii) &
                          - 1/(1D0-km**2)* cos(km*ii*pi/im)
       endif
       g_X_Weight(ii) = 2D0/im * g_X_Weight(ii) * xl/2
    enddo
    g_X_Weight(0)  = g_X_Weight(0) / 2
    g_X_Weight(im) = g_X_Weight(im) / 2

    call MessageNotify('M','at_initial','at_module (2009/07/31) is initialized')

  end subroutine at_Initial

! ---- 逆変換 ----
  function ag_at(at_data)
    !
    ! チェビシェフデータから格子データへ変換する(2 次元配列用).
    !
    real(8), dimension(:,:), intent(in)      :: at_data
    !(in) チェビシェフデータ

    real(8), dimension(size(at_data,1),0:im) :: ag_at
    !(out) 格子点データ

    real(8), dimension(size(at_data,1)*im) :: y
    ! 作業用配列
    integer :: m

    m = size(at_data,1)
    if ( size(at_data,2)-1 < km ) then
       call MessageNotify('E','ag_at', &
            'The Chebyshev dimension of input data too small.')
    elseif ( size(at_data,2)-1 > km ) then
       call MessageNotify('W','ag_at', &
            'The Chebyshev dimension of input data too large.')
    endif

    ag_at = 0.0D0
    ag_at(:,0:km)=at_data
    call fttctb(m,im,ag_at,y,it,t)

  end function ag_at

  function g_t(t_data)
    !
    ! チェビシェフデータから格子データへ変換する(1 次元配列用).
    !
    real(8), dimension(:), intent(in)  :: t_data
    !(in) チェビシェフデータ

    real(8), dimension(0:im)           :: g_t
    !(out) 格子点データ

    real(8), dimension(1,size(t_data)) :: t_work
    ! 作業用配列
    real(8), dimension(1,0:im)         :: g_work
    ! 作業用配列

    t_work(1,:) = t_data
    g_work = ag_at(t_work)
    g_t = g_work(1,:)

  end function g_t


! ---- 正変換 ----
  function at_ag(ag_data)
    !
    ! 格子データからチェビシェフデータへ変換する(2 次元配列用).
    !
    real(8), dimension(:,:), intent(in)      :: ag_data
    !(in) 格子点データ

    real(8), dimension(size(ag_data,1),0:km) :: at_ag
    !(out) チェビシェフデータ

    real(8), dimension(size(ag_data,1)*im)   :: y
    real(8), dimension(size(ag_data,1),0:im) :: ag_work
    integer :: m

    m = size(ag_data,1)
    if ( size(ag_data,2)-1 < im ) then
       call MessageNotify('E','at_ag', &
            'The Grid points of input data too small.')
    elseif ( size(ag_data,2)-1 > im ) then
       call MessageNotify('W','at_ag', &
            'The Grid points of input data too large.')
    endif
    ag_work = ag_data

    call fttctf(m,im,ag_work,y,it,t)
    at_ag = ag_work(:,0:km)

  end function at_ag

  function t_g(g_data)  ! 台形格子 -> スペクトル
    !
    ! 格子データからチェビシェフデータへ変換する(1 次元配列用).
    !
    real(8), dimension(:), intent(in)     :: g_data
    !(in) 格子点データ

    real(8), dimension(0:km)              :: t_g
    !(out) チェビシェフデータ

    real(8), dimension(1,size(g_data)) :: ag_work
    real(8), dimension(1,0:km)         :: at_work

    ag_work(1,:) = g_data
    at_work = at_ag(ag_work)
    t_g = at_work(1,:)

  end function t_g

! ---- 微分計算 ----
  function at_Dx_at(at_data)
    !
    ! 入力チェビシェフデータに X 微分を作用する(2 次元配列用).
    !
    ! チェビシェフデータの X 微分とは, 対応する格子点データに X 微分を
    ! 作用させたデータのチェビシェフ変換のことである.
    !
    !
    real(8), dimension(:,0:), intent(in)                    :: at_data
    !(in) 入力チェビシェフデータ

    real(8), dimension(size(at_data,1),0:size(at_data,2)-1) :: at_Dx_at
    !(out) チェビシェフデータの X 微分

    integer :: m, k
    integer :: nm, kmax

    nm=size(at_data,1)
    kmax=size(at_data,2)-1
    if ( kmax  < km ) then
       call MessageNotify('W','at_Dx_at', &
            'The Chebyshev dimension of input data too small.')
    elseif ( kmax > km ) then
       call MessageNotify('E','at_Dx_at', &
            'The Chebyshev dimension of input data too large.')
    endif

    if ( kmax == im ) then
!       !$omp parallel do
       do m=1,nm
          at_Dx_at(m,kmax)   = 0.
          at_Dx_at(m,kmax-1) = 2 * km * at_data(m,kmax) /2
       enddo
    else
       do m=1,nm
          at_Dx_at(m,kmax)   = 0.
          at_Dx_at(m,kmax-1) = 2 * km * at_data(m,kmax)
          ! スタートはグリッド対応最大波数未満. Factor 1/2 不要
       enddo
    endif

    do k=kmax-2,0,-1
  !     !$omp parallel do
       do m=1,nm
          at_Dx_at(m,k) = at_Dx_at(m,k+2) + 2*(k+1)*at_data(m,k+1)
       enddo
    enddo

 !   !$omp parallel do private(m)
    do k=0,kmax
       do m=1,nm
          at_Dx_at(m,k) = 2/xl * at_Dx_at(m,k)
       enddo
    enddo

  end function at_Dx_at

  function t_Dx_t(t_data)
    !
    ! 入力チェビシェフデータに X 微分を作用する(1 次元配列用).
    !
    ! チェビシェフデータの X 微分とは, 対応する格子点データに X 微分を
    ! 作用させたデータのチェビシェフ変換のことである.
    !
    !
    real(8), dimension(:), intent(in)   :: t_data
    !(in) 入力チェビシェフデータ

    real(8), dimension(size(t_data))    :: t_Dx_t
    !(out) チェビシェフデータの X 微分

    real(8), dimension(1,size(t_data))  :: at_work
    ! 作業用配列

    at_work(1,:) = t_data
    at_work = at_Dx_at(at_work)
    t_Dx_t = at_work(1,:)

  end function t_Dx_t

! ---- 補間計算 ----
  function Interpolate_t(t_data,xval)
    real(8), dimension(0:), intent(in)  :: t_data
    !(in) 入力チェビシェフデータ

    real(8), intent(in)                 :: xval
    ! 補間する点の座標

    real(8)                             :: Interpolate_t
    ! 補間した結果の値

    integer :: kmax
    ! 入力配列の最大次数

    real(8) :: y2, y1, y0, x
    ! Crenshow's reccurence formula 計算用変数

    integer :: k
    ! DO 文変数

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
    !(in) 入力チェビシェフデータ

    real(8), intent(in)                  :: xval
    ! 補間する点の座標

    real(8), dimension(size(at_data,1))  :: a_Interpolate_at
    ! 補間した結果の値

    integer :: kmax
    ! 入力配列の最大次数

    real(8), dimension(size(at_data,1))  :: y2, y1, y0
    real(8)                              ::  x
    ! Crenshow's reccurence formula 計算用変数

    integer :: k
    ! DO 文変数

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

!---- Dirichlet 型境界条件(タウ法) ----

  subroutine at_BoundariesTau_DD_2d(at_data,values)
    !
    ! Dirichlet 型境界条件の適用(タウ法, 2 次元配列用)
    ! 両境界での値を与える.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data ! データ(m,0:km)
    !(inout) 境界条件を適用するチェビシェフデータ(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) 境界値(m,2)

    real(8), dimension(:,:), allocatable  :: alu
    integer, dimension(:), allocatable    :: kp
    real(8), dimension(0:km,0:km)         :: tt_data
    real(8), dimension(0:km,0:im)         :: tg_data
    real(8), dimension(size(at_data,1))    :: value1, value2           ! 境界値

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
    ! Dirichlet 型境界条件の適用(タウ法, 1 次元配列用)
    ! 両境界での値を与える.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) 境界条件を適用するチェビシェフデータ(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) 境界値

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! 境界値

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesTau_DD_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesTau_DD_1d


!---- Dirichlet/Neumann 型境界条件(タウ法) ----

  subroutine at_BoundariesTau_DN_2d(at_data,values)
    !
    ! Dirichlet/Neumann 型境界条件の適用(タウ法, 2 次元配列用)
    ! i=0 で値, i=im で勾配の値を与える.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) 境界条件を適用するチェビシェフデータ(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) 境界値(m,2)

    real(8), dimension(:,:), allocatable  :: alu
    integer, dimension(:), allocatable    :: kp
    real(8), dimension(0:km,0:km)         :: tt_data
    real(8), dimension(0:km,0:im)         :: tg_data
    real(8), dimension(size(at_data,1))    :: value1, value2           ! 境界値

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
    ! Dirichlet/Neumann 型境界条件の適用(タウ法, 1 次元配列用)
    ! i=0 で値, i=im で勾配の値を与える.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) 境界条件を適用するチェビシェフデータ(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) 境界値

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! 境界値

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesTau_DN_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesTau_DN_1d

!---- Neumann/Dirichlet 型境界条件(タウ法) ----

  subroutine at_BoundariesTau_ND_2d(at_data,values)
    !
    ! Neumann/Dirichlet 型境界条件の適用(タウ法, 2 次元配列用)
    ! i=0 で勾配の値, i=im で値を与える.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) 境界条件を適用するチェビシェフデータ(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) 境界値(m,2)

    real(8), dimension(:,:), allocatable  :: alu
    integer, dimension(:), allocatable    :: kp
    real(8), dimension(0:km,0:km)         :: tt_data
    real(8), dimension(0:km,0:im)         :: tg_data
    real(8), dimension(size(at_data,1))    :: value1, value2           ! 境界値

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
    ! Neumann/Dirichlet 型境界条件の適用(タウ法, 1 次元配列用)
    ! i=0 で勾配の値, i=im で値を与える.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) 境界条件を適用するチェビシェフデータ(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) 境界値

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! 境界値

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesTau_ND_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesTau_ND_1d


!---- Neumann 型境界条件(タウ法) ----

  subroutine at_BoundariesTau_NN_2d(at_data,values)
    !
    ! Neumann/Dirichlet 型境界条件の適用(タウ法, 2 次元配列用)
    ! 両境界で勾配の値を与える.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) 境界条件を適用するチェビシェフデータ(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) 境界値(m,2)

    real(8), dimension(:,:), allocatable  :: alu
    integer, dimension(:), allocatable    :: kp
    real(8), dimension(0:km,0:km)         :: tt_data
    real(8), dimension(0:km,0:im)         :: tg_data
    real(8), dimension(size(at_data,1))    :: value1, value2           ! 境界値

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
    ! Neumann/Dirichlet 型境界条件の適用(タウ法, 1 次元配列用)
    ! 両境界で勾配の値を与える.
    !
    ! このサブルーチンを直接使うことを勧めない.
    ! 共通インターフェース at_Boundaries_NN を用いること.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) 境界条件を適用するチェビシェフデータ(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) 境界値

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! 境界値

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesTau_NN_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesTau_NN_1d

  !--------------- 積分計算 -----------------
    function a_Int_ag(ag)
      !
      ! 1 次元格子点データが並んだ 2 次元配列の積分
      !
      real(8), dimension(:,0:), intent(in)     :: ag
      !(in)入力格子点データ

      real(8), dimension(size(ag,1))           :: a_Int_ag
      !(out) 積分したデータ
      integer :: i

      if ( size(ag,2) < im+1 ) then
         call MessageNotify('E','ae_ag', &
              'The Grid points of input data too small.')
      elseif ( size(ag,2) > im+1 ) then
         call MessageNotify('W','ae_ag', &
              'The Grid points of input data too large.')
      endif

      a_Int_ag = 0.0d0

 !     !$omp parallel do reduction(+:a_Int_ag)
      do i=0,im
         a_Int_ag(:) = a_Int_ag(:) + ag(:,i)*g_X_Weight(i)
      enddo
    end function a_Int_ag

    function Int_g(g)
      !
      ! 1 次元格子点データの積分および平均.
      !
      real(8), dimension(0:im), intent(in)   :: g
      !(in) 格子点データ

      real(8)                                :: Int_g
      !(out) 積分値

      Int_g = sum(g*g_X_Weight)
    end function Int_g

    function a_Avr_ag(ag)
      !
      ! 1 次元格子点データが並んだ 2 次元配列の平均
      !
      real(8), dimension(:,0:), intent(in)   :: ag
      !(in)入力格子点データ

      real(8), dimension(size(ag,1))         :: a_Avr_ag
      !(out) 平均したデータ

      a_Avr_ag = a_Int_ag(ag)/sum(g_X_Weight)
    end function a_Avr_ag

    function Avr_g(g)
      !
      ! 1 次元格子点データの平均
      !
      real(8), dimension(0:im), intent(in)   :: g
      !(in) 格子点データ

      real(8)                                :: Avr_g
      !(out) 積分値

      Avr_g = Int_g(g)/sum(g_X_Weight)
    end function Avr_g



!---- Dirichlet 型境界条件(実空間での評価) ----

  subroutine at_BoundariesGrid_DD_2d(at_data,values)
    !
    ! Dirichlet 型境界条件の適用(実空間での評価, 2 次元配列用)
    ! 両境界での値を与える.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) 境界条件を適用するチェビシェフデータ(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) 境界値(m,2)

    real(8), dimension(:,:), allocatable     :: alu
    integer, dimension(:), allocatable       :: kp
    real(8), dimension(size(at_data,1),0:im) :: ag_data
    real(8), dimension(0:km,0:km)            :: tt_data
    real(8), dimension(0:km,0:im)            :: tg_data
    real(8), dimension(size(at_data,1))      :: value1, value2  ! 境界値

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
    ! Dirichlet 型境界条件の適用(実空間での評価, 1 次元配列用)
    ! 両境界での値を与える.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) 境界条件を適用するチェビシェフデータ(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) 境界値

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! 境界値

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesGrid_DD_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesGrid_DD_1d

!---- Dirichlet/Neumann 型境界条件(実空間での評価) ----

  subroutine at_BoundariesGrid_DN_2d(at_data,values)
    !
    ! Dirichlet/Neumann 型境界条件の適用(実空間での評価, 2 次元配列用)
    ! i=0 で値, i=im で勾配の値を与える.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) 境界条件を適用するチェビシェフデータ(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) 境界値(m,2)

    real(8), dimension(:,:), allocatable     :: alu
    integer, dimension(:), allocatable       :: kp
    real(8), dimension(size(at_data,1),0:im) :: ag_data
    real(8), dimension(0:km,0:km)            :: tt_data
    real(8), dimension(0:km,0:im)            :: tg_data
    real(8), dimension(size(at_data,1))      :: value1, value2  ! 境界値

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
    ! Dirichlet/Neumann 型境界条件の適用(実空間での評価, 1 次元配列用)
    ! i=0 で値, i=im で勾配の値を与える.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) 境界条件を適用するチェビシェフデータ(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) 境界値

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! 境界値

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesGrid_DN_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesGrid_DN_1d

!---- Neumann/Dirichlet 型境界条件(実空間での評価) ----

  subroutine at_BoundariesGrid_ND_2d(at_data,values)
    !
    ! Neumann/Dirichlet 型境界条件の適用(実空間での評価, 2 次元配列用)
    ! i=0 で勾配の値, i=im で値を与える.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) 境界条件を適用するチェビシェフデータ(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) 境界値(m,2)

    real(8), dimension(:,:), allocatable     :: alu
    integer, dimension(:), allocatable       :: kp
    real(8), dimension(size(at_data,1),0:im) :: ag_data
    real(8), dimension(0:km,0:km)            :: tt_data
    real(8), dimension(0:km,0:im)            :: tg_data
    real(8), dimension(size(at_data,1))      :: value1, value2  ! 境界値

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
    ! Neumann 型境界条件の適用(実空間での評価, 1 次元配列用)
    ! i=0 で勾配の値, i=im で値を与える.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) 境界条件を適用するチェビシェフデータ(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) 境界値

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! 境界値

    if (.not. present(values)) then
       vwork(1,1)=0 ; vwork(1,2)=0
    else
       vwork(1,:) = values
    endif

    at_work(1,:)=t_data
    call at_BoundariesGrid_ND_2d(at_work,vwork)
    t_data=at_work(1,:)

  end subroutine at_BoundariesGrid_ND_1d

!---- Neumann 型境界条件 ----

  subroutine at_BoundariesGrid_NN_2d(at_data,values)
    !
    ! Neumann 型境界条件の適用(実空間での評価, 2 次元配列用)
    ! 両境界で勾配の値を与える.
    !
    real(8), dimension(:,0:),intent(inout)         :: at_data
    !(inout) 境界条件を適用するチェビシェフデータ(m,0:km)

    real(8), dimension(:,:), intent(in), optional  :: values
    !(in) 境界値(m,2)

    real(8), dimension(:,:), allocatable     :: alu
    integer, dimension(:), allocatable       :: kp
    real(8), dimension(size(at_data,1),0:im) :: ag_data
    real(8), dimension(0:km,0:km)            :: tt_data
    real(8), dimension(0:km,0:im)            :: tg_data
    real(8), dimension(size(at_data,1))      :: value1, value2  ! 境界値

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
    ! Neumann 型境界条件の適用(実空間での評価, 1 次元配列用)
    ! 両境界で勾配の値を与える.
    !
    real(8), dimension(0:km),intent(inout)       :: t_data
    !(inout) 境界条件を適用するチェビシェフデータ(0:km)

    real(8), dimension(2), intent(in), optional  :: values
    !(in) 境界値

    real(8), dimension(1,0:km)                   :: at_work
    real(8), dimension(1,2)                      :: vwork           ! 境界値

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
