=begin
#
# * Yuta Kawai
#   * 過去の更新履歴
#      * 2018/11/30(河合佑太) INSTALL ドキュメント作成
=end

=begin JA
= DCPOM インストールガイド
=end JA

=begin EN
= DCPOM Installation Guide
=end EN

=begin JA
== 動作環境
DCPOM は以下の環境での動作を確認しています.
=end JA

=begin EN
== Operation Environment
DCPOM is operated by following environments.
=end EN

=begin

* OS
  * SUSE Linux Enterprise Server 11 SP3
  * Debian GNU/Linux ver. 8.11
  * MacOS High Sierra ver. 10.13.6
* Fortran Compiler
  * GFortran ver. 8.2.0 
    * with LAPACK  ver. 3.5.0
    * with OpenMPI ver. 3.1.3
  * GFortran ver. 4.9.2 
    * with LAPACK  ver. 3.8.0
    * with OpenMPI ver. 1.6.5    
  * Intel &reg; Compilers ver. 16.0.1
    * with Intel &reg MKL 
    * with Intel &reg MPI Libray ver. 5.1.2

=end

=begin JA
== 必要なソフトウェア

DCPOM を利用するためには, 以下のソフトウェアを
事前にインストールしておく必要があります.

* ((<netCDF|URL:http://www.gfd-dennou.org/library/netcdf>)),  >= ver. 3.6
* ((<ISPACK|URL:http://www.gfd-dennou.org/library/ispack>)), ver. 1.0.4
* ((<Gtool5|URL:https://www.gfd-dennou.org/library/gtool/gtool5.htm.ja>)), >= ver. 20101218-1
* ((<SPMODEL|URL:http://www.gfd-dennou.org/library/spmodel>)), >= ver. 0.8.0
* (MPI library)

=end

=begin EN
== Software Requirements

The following software needs to use DCPOM
* ((<netCDF|URL:http://www.gfd-dennou.org/library/netcdf>)),  >= ver. 3.6
* ((<ISPACK|URL:http://www.gfd-dennou.org/library/ispack/index.htm.en>)), ver. 1.0.4
* ((<Gtool5|URL:http://www.gfd-dennou.org/library/gtool/gtool5.htm.en>)), >= ver. 20101218-1
* ((<SPMODEL|URL:http://www.gfd-dennou.org/library/spmodel/index.htm.en>)), ver. 0.8.0
* (MPI library)

=end EN

=begin JA
== ビルドの手引き
=end JA

=begin EN
== How to build
=end EN

=begin JA
=== TGZ パッケージの展開

適当な作業ディレクトリでソースアーカイブを展開します.
ソースは dcpom-((|バージョン|)) というディレクトリに展開されます.

        $ tar xvzf dcpom_current.tgz

または

        $ zcat dcpom_current.tar.gz | tar -xvf -

=end JA

=begin EN
=== Extract TGZ Package

Make an empty directory, and extract archive.
A directory `dcpom-((|version|))'
created at the current working directory.

        $ tar xvzf dcpom_current.tgz

or

        $ zcat dcpom_current.tar.gz | tar -xvf -

=end EN

=begin JA
=== Fortran コンパイラの指定

環境変数 ((* FC *)) に使用する Fortran コンパイラを指定してください.
以下は, 利用するコンパイラが gfortran の場合です.

* sh, bash の場合

        $ FC=gfortran ; export FC

* csh, tcsh の場合

        $ setenv FC gfortran

最適化やデバッグのためのオプションは環境変数 ((* FCFLAGS *))
に設定してください. 以下の例は gfortran を使用する場合の最適化と
スレッド並列化のためのオプションです.

* sh, bash の場合

        $ FCFLAGS="-O3 -fopenmp" ; export FCFLAGS

* csh, tcsh の場合

        $ setenv FCFLAGS "-O3 -fopenmp"


=end JA
=begin EN
=== Specify Fortran Compiler

Specify Fortran compiler to environment variable ((* FC *)).
For example, if you use "gfortran",

* sh, bash

        $ FC=gfortran ; export FC

* csh, tcsh

        $ setenv FC gfortran

Specify Fortran compiler options for optimization and debug to
environment variable ((* FCFLAGS *)).
For example, if you set options for optimization and
thread parallelization to gfortran,

* sh, bash

        $ FCFLAGS="-O3 -fopenmp" ; export FCFLAGS

* csh, tcsh

        $ setenv FCFLAGS "-O3 -fopenmp"


=end EN


=begin JA
=== configure の実行

展開されたディレクトリに移動し, (({ ./configure }))を実行します. 
その際, --with-*= で使用するライブラリのパスを指定します. 

* (({ --with-netcdf= })) には netCDF ライブラリのパスを指定します.
  (以下の例は /usr/local/netcdf/lib/libnetcdf.a にライブラリがある
  場合のものです).
* (({ --with-ispack= })) には ispack ライブラリのパスを指定します.
  (以下の例は /usr/local/ispack/lib/libisp.a にライブラリがある
  場合のものです).
* (({ --with-gtool5= })) には gtool5 ライブラリのパスを指定します.
  (以下の例は /usr/local/gtool5/lib/libgtool5.a にライブラリがある
  場合のものです).
* (({ --with-spml= })) には SPOMODEL ライブラリのパスを指定します.
  (以下の例は /usr/local/ispack/lib/libspml-omp.a にライブラリがある
  場合のものです).
    
(({ --with-netcdff= })) も指定する必要があるかもしれません(詳しくは下記のオプションの詳細を参照).
以下のコマンドを実行すると, Makefile が作成されます. 

        $ ./configure --with-netcdf=/usr/local/netcdf/lib/libnetcdf.a \
                      --with-ispack=/usr/local/ispack/lib/libisp.a    \
                      --with-gtool5=/usr/local/gtool5/lib/libgtool5.a \
                      --with-spml=/usr/local/spml/lib/libspml-omp.a

DCPOM を MPI 用にビルドする場合には,((<MPI を使用した計算用にビルドする場合>)) を参照してください.
また, 軸対称計算を行いたい場合は「--with-dogcm-mode=AXISYM」を,  
緯度方向のスペクトル変換に sjpack を使用したい場合は「--enable-sjapck」を追加してください. 
以下のように (({ --help })) オプションをつけることで, その他の指定可能なオプションリストが表示されます.

        $ ./configure --help

主なオプションに関しての説明です.

:(({--with-netcdf=}))((|ARG|))
  ((|ARG|)) に ((<ビルドに必要な netCDF ライブラリ>))
  を指定します. 必ず明示的に指定する必要があります.

:(({--with-netcdff=}))((|ARG|))
  netCDF ライブラリが共有ライブラリである場合, C 用ライブラリと
  Fortran 用ライブラリとに分かれてビルドされている場合があります.
  その際は, 上記オプションに C 用ライブラリを指定し, 本オプションの
  ((|ARG|)) に ((<Fortran 用ライブラリ>)) を指定します.

:(({--with-netcdf-include=}))((|ARG|))
  必要ならば netCDF ライブラリの Fortran 用ヘッダ(netcdf.inc) を指定します.

:(({--with-gtool5=}))((|ARG|))
  ((|ARG|)) に ((<ビルドに必要な Gtool5 ライブラリ>))
  を指定します. 必ず明示的に指定する必要があります.

:(({--with-ispack=}))((|ARG|))
  ((|ARG|)) に ((<ビルドに必要な ISPACK ライブラリ>))
  を指定します. 必ず明示的に指定する必要があります.

:(({--with-spml=}))((|ARG|))
  ((|ARG|)) に ((<ビルドに必要な SPMODEL ライブラリ>))
  を指定します. 必ず明示的に指定する必要があります.

:(({--with-lapack=}))((|ARG|))
  ((|ARG|)) に ((<連立一次方程式を解くために使用するライブラリ>))
  を指定します.

:(({--prefix=}))((|ARG|))
  ((|ARG|)) にライブラリやモジュール, 実行ファイルのインストール先の
  ディレクトリのプレフィックスを指定します.
  デフォルトは (({ /usr/local/spml })) です.

:(({--host=}))((|ARG|))
  クロスコンパイルを行う場合には, パッケージが実行されるシステムタイプ名
  を ((|ARG|)) に指定します.

:(({--libdir=}))((|ARG|))
  ((|ARG|)) にライブラリのインストール先のディレクトリを指定します.
  デフォルトは (({ /usr/local/dcpom/lib })) です.

:(({--includedir=}))((|ARG|))
  ((|ARG|)) にモジュール情報ファイルのインストール先のディレクトリ
  を指定します. デフォルトは (({ /usr/local/dcpom/include })) です.

:(({--bindir=}))((|ARG|))
  ((|ARG|)) に実行ファイルのインストール先のディレクトリを指定します.
  デフォルトは (({ /usr/local/dcpom/bin })) です.

:(({--docdir=}))((|ARG|))
  ((|ARG|)) にドキュメント/マニュアルのインストール先のディレクトリを指定します.
  デフォルトは (({ /usr/local/dcpom/doc })) です.

=end JA
=begin EN

=== Execute `configure'

Move created directroy, and excute `(({ ./configure }))'.
You should set necessary library location, netCDF, ISPACK, Gtool5 and SPMODEL.
If your path of netCDF library is `/usr/local/netcdf/lib/libnetcdf.a',
ISPACK library is `/usr/local/ispack/lib/libisp.a',
Gtool5 library is `/usr/local/gtool5/lib/libgtool5.a', and
SPMODEL library is `/usr/local/spmodel/lib/libspml-omp.a',
you should set options as follow.
Then a `Makefile' will be created at
the current working directory.
If the netCDF library is a shared library, (({ --with-netcdff= }))
option may be needed.
See details of options as follows.

        $ ./configure --with-netcdf=/usr/local/netcdf/lib/libnetcdf.a \
                --with-ispack=/usr/local/ispack/lib/libisp.a    \
                --with-gtool5=/usr/local/gtool5/lib/libgtool5.a \
                --with-spml=/usr/local/spml/lib/libspml-omp.a

If DCPOM is built for MPI, see ((<How to build for MPI>)).
For axisymmetric calculations, add a option, "--with-dogcm-mode=AXISYM". 
If users want to use sjapck for spectral transformation, add a option, "--enable-sjapck". 

In the case of setting (({ --help })) option as follow, Available
options are showed.

        $ ./configure --help

Descriptions about principal options are listed below.

:(({--with-netcdf=}))((|ARG|))
  Specify ((<netCDF library needed for build>)) to ((|ARG|)).
  You must specify explicitly.

:(({--with-netcdff=}))((|ARG|))
  If the netCDF library is a shared library, the library may be divided
  C library from Fortran library. In the case, specify the C library
  to above option, and specify
  ((<the Fortran library>)) to ((|ARG|)) in this option.

:(({--with-netcdf-include=}))((|ARG|))
  Set location of netCDF header file for fortran(netcdf.inc), if you need.

:(({--with-ispack=}))((|ARG|))
  Specify ((<ISPACK library needed for build>)) to ((|ARG|)).
  You must specify explicitly.
    
:(({--with-gtool5=}))((|ARG|))
  Specify ((<Gtool5 library needed for build>)) to ((|ARG|)).
  You must specify explicitly.

:(({--with-spml=}))((|ARG|))
  Specify ((<SPMODEL library needed for build>)) to ((|ARG|)).
  You must specify explicitly.
    
:(({--with-lapack=}))((|ARG|))
  Specify ((<library using solve liner equations >)) to ((|ARG|)),
  if you need.

:(({--prefix=}))((|ARG|))
  Specify prefix to ((|ARG|)).
  Default value is (({ /usr/local/dcpom })).

:(({--host=}))((|ARG|))
  When cross-compiling, specify
  the type of system on which the package will run to
  ((|ARG|)).

:(({--libdir=}))((|ARG|))
  Specify directory to which the library is installed to ((|ARG|)).
  Default value is (({ /usr/local/dcpom/lib })).

:(({--includedir=}))((|ARG|))
  Specify directory to which the module is installed to ((|ARG|)).
  Default value is (({ /usr/local/dcpom/include })).

:(({--bindir=}))((|ARG|))
  Specify directory to which the executable file is installed to ((|ARG|)).
  Default value is (({ /usr/local/dcpom/bin })).

:(({--with-docdir=}))((|ARG|))
  Specify directory to which the documentation file is installed to ((|ARG|)).
  Default value is (({ /usr/local/dcpom/doc })).

=end EN

=begin JA
=== MPI を使用した計算用にビルドする場合

MPI 並列化した DPCOM を使用したい場合は, MPI ライブラリのインストールに加えて, 
MPI サポートを有効にして Gtool5 や SPMODEL を インストールする必要があります. 
詳細については
((<Gtool5 のインストールドキュメント|URL:http://www.gfd-dennou.org/library/gtool/gtool5/gtool5_current/INSTALL.htm>)), 
((<SPMODEL のインストールドキュメント|URL:http://www.gfd-dennou.org/library/spmodel/spml_current/doc/INSTALL.htm.ja>)), 
を参照して下さい.

また, DCPOM の((<configure の実行>))を行うときには, 
環境変数 FC に mpif90 などの MPI 用のコンパイルコマンドを指定する必要があります. 

=end JA

=begin EN
=== How to build for calculations with MPI

If users want to use parallelized DCPOM with MPI, in addition to installation of MPI library, 
installation of Gtool5 and SPMODEL should be performed so that MPI is activated.  
For the detail, see
((<"Gtool5 Installation Guide"|URL:http://www.gfd-dennou.org/library/gtool/gtool5/gtool5_current/INSTALL.htm.en>)).
and ((<"SPMODEL Installation Guide"|URL:http://www.gfd-dennou.org/library/spmodel/spml_current/doc/INSTALL.htm.en>)).

When executing "configure" for DCPOM, specify a compile commend like as "mpif90"
to environment variable ((* FC *)). 

=end EN

=begin JA
=== ソースコードのコンパイル

configure を実行すると Makefile が更新されます.
make を実行してコンパイルを行なって下さい.

  $ make

「libsrc, tool, model/globalSWModel_DG, model/globalModel_FVM, model/dogcm」ディレクトリの順番に, 
ソースコードがコンパイルされます. 

=end JA
=begin EN
=== Compile source code

When ./configure is executed, Makefile is updated
displayed as follows. Execute `make' as follows:

  $ make

Source codes will be compiled in order of "libsrc, tool, model/globalSWModel_DG, 
model/globalModel_FVM, model/dogcm" directories. 

=end EN

=begin JA
=== ドキュメントの生成

モデルのリファレンスやユーザガイドのコンパイルは, トップディレクトリにおいて, 以下のコマンドを実行してください.
((<DCPOMの TGZ パッケージ|URL:http://www.gfd-dennou.org/library/dcpom/dcpom_latest.tgz>))
から入手する場合には既に生成済みです. 
なお, ドキュメントの生成には, ((<rdtool|URL://https://uwabami.github.io//software/rdtool/index.html>))等の
幾つかのソフトウェアが必要です. 

        $ make update-doc
=end JA
=begin EN
=== Generate documentations

To generate documentations, execute following command in current directory.
If you get from
((<DCPOM TGZ package|URL:http://www.gfd-dennou.org/library/dcpom/dcpom/dcpom_latest.tgz>)),
documentations are already generated.
Note that some softwares such as ((<rdtool|URL://https://uwabami.github.io//software/rdtool/index.html>)) 
are necessary to generate documentations. 

        $ make udpate-doc
=end EN

=begin JA
== プログラム実行の手順

DCPOM のビルドが成功した場合は, 例えば, model/dogcm ディレクトリの下に,
dogcm (DCPOM に含まれる海洋海氷モデルの一つ)や実験結果の解析用プログラムが,  
((*dogcm_axisym*)) (軸対象設定の場合), ((*ocndiag*)) として作成されます. 

これらの実行の手順については, 
((<DCPOM ユーザーガイド|URL:users-guide/dcpom_users-guide.pdf>)) を参照してください. 

=end JA
=begin EN
== Execute programs

If DCPOM is built successfully, for example, two binaries for a ocean-sea ice model in DCPOM, dogcm, and 
a program to analysis the simulation results are generated under "model/dogcm", 
as ((*dogcm_axisym*)) (in case of axisymmetric setting) and ((*ocndiag*)), respectively. 

See ((<DCPOM users guide|URL:users-guide/dcpom_users-guide.pdf>)) for information about the execution procedures. 

=end EN

