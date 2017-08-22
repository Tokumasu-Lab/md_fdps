# 開発環境構築メモ
川井喜与人 2015/11/4  初版作成  
川井喜与人 2015/12/28 追記  
川井喜与人 2016/1/6   補足:MPIプログラムのデバッグ をGDBの使い方全般として独立 ->  "README_GDB_usage.txt"  
川井喜与人 2016/1/31  GCCのバージョンを5.3.0に変更． 複数のGCCの切り替え方について追記  
川井喜与人 2016/3/3   GCCのバージョンとパスの確認方法について追記  
川井喜与人 2016/3/9   OpenMPIのバージョン確認方法について追記  
川井喜与人 2016/8/16   GCC,OpenMPIのコンパイルオプションを修正(GCC-6.1.0, OpenMPI-2.0.0で再検証)  
川井喜与人 2016/12/14  Windows Explorer からソースファイル内の全文検索(フォルダごと)を有効化する手順を追記  
川井喜与人 2017/4/7  メモ全体を markdown 形式で再作成．全体的に加筆，修正  
川井喜与人 2017/4/8  おまけ2 のスクリプトを使用すると sftp 接続ができなくなる問題を修正  
川井喜与人 2017/8/23 公開用に内容を修正


### この文書の内容
- [現在(2017/4/7)のIFSスパコンの環境と，ユーザーが新たに構築する環境について](#IFSスパコンの環境)
- [インストールの準備](#インストールの準備)
- [GCC のコンパイルとインストール](#GCCのコンパイルとインストール)
- [OpenMPI のコンパイルとインストール](#OpenMPIのコンパイルとインストール)
- [FFTw のコンパイルとインストール](#FFTwのコンパイルとインストール)
- [おまけ1：　FDPSのParticleMesh機能拡張のコンパイル](#FDPSのparticle_mesh拡張機能のコンパイル)
- [おまけ2：　複数のGCCバージョンの切り替え](#複数のgccバージョンの切り替え)

<a id="IFSスパコンの環境"></a>
<a href="#IFSスパコンの環境"></a>
## 現在(2017/4/7)のIFSスパコンの環境と，ユーザーが新たに構築する環境について
FDPS (Framework for Developing Particle Simulator) を使用して分子動力学プログラムを作成するため，FDPSを利用したプログラムをコンパイルするために必要な環境を構築する．

FDPS ver3.0 向け  
参考：
- https://github.com/FDPS/FDPS
- http://www.jmlab.jp/?p=530

#### 開発元の動作確認環境
C++ コンパイラ (gcc 4.4.5以降，あるいはKコンパイラ 1.2.0)  
MPI version 1.3環境　(OpenMPI 1.8.1で動作確認)  
FFTw 3.3以降  

現在(2017/4/7)のIFSスパコン(Phase-2: afisfep2)の環境  
- シェル bash
- gcc 4.3.4  
MPI環境:  
  - OpenMPI 1.5.3  
  - SGI MPT 2.13  
  - Intel MPI 4.1.1.036  

- MKL 11.0.5.192  
(intelの数学ライブラリ．FFTW互換のインターフェイスはあるが完全な互換性は保障されない)

GCCが古すぎるのもあるが，C++11以降の規格が使用できると便利なこともあり，開発環境として以下の環境をユーザのローカルフォルダにインストールした．
- gcc 6.3.0
- OpenMPI 2.0.1
- FFTW 3.3.6


<a id="インストールの準備"></a>
<a href="#インストールの準備"></a>
## インストールの準備
適当な作業フォルダ(ここでは `$HOME/install/` )に上記各ソフトウェアの圧縮ソースファイルをダウンロードし，展開しておく．  
下記で例示しているディレクトリ( `/home/hogehoge` )は作業者のホームディレクトリである．各自 `$ pwd` コマンド等で確認すること．  

```
$ pwd
/home/hogehoge/install

$ ls
fftw-3.3.6-pl1.tar.gz
gcc-6.3.0.tar.bz2
openmpi-2.0.1.tar.bz2

$ tar -xvf gcc-6.3.0.tar.bz2
$ tar -xvf openmpi-2.0.1.tar.bz2
$ tar -zxvf fftw-3.3.6-pl1.tar.gz
```

既存環境(ここではgcc-4.3)との切り替えのため，コマンドをシンボリックリンクとして作成する．  
パスを登録するディレクトリ(ここでは `$HOME/local_path` )を作成

```
$ cd
$ mkdir local_path
$ cd local_path
```

既存環境のコンパイラのシンボリックリンクの作成

```
$ pwd
/home/hogehoge/local_path

$ ln -s /usr/bin/gcc gcc-4.3
$ ln -s /usr/bin/g++ g++-4.3
$ ln -s /usr/bin/gfortran gfortran-4.3
```

シンボリックリンクのパスの登録  
.bashrc の最後に以下を追加．　あるいは[おまけ2](#複数のgccバージョンの切り替え)のスクリプトを利用する．
```
export PATH=$HOME/local_path:$PATH
```
.bashrcを変更後，__ログインしなおして__ 登録したシンボリックリンク(例えば `gcc-4.3` )が有効であることを確認する．
```
$ gcc-4.3 -v
Using built-in specs.
Target: x86_64-suse-linux
Configured with: ../configure --prefix=/usr --infodir=/usr/share/info --mandir=/usr/share/man --libdir=/usr/lib64 --libexecdir=/usr/lib64 --enable-languages=c,c++,objc,fortran,obj-c++,java,ada --enable-checking=release --with-gxx-include-dir=/usr/include/c++/4.3 --enable-ssp --disable-libssp --with-bugurl=http://bugs.opensuse.org/ --with-pkgversion='SUSE Linux' --disable-libgcj --disable-libmudflap --with-slibdir=/lib64 --with-system-zlib --enable-__cxa_atexit --enable-libstdcxx-allocator=new --disable-libstdcxx-pch --enable-version-specific-runtime-libs --program-suffix=-4.3 --enable-linux-futex --without-system-libunwind --with-cpu=generic --build=x86_64-suse-linux
Thread model: posix
```

コンパイル作業用のディレクトリを作成  

```
$ cd
$ cd install

$ mkdir build
$ cd build

$ pwd
/home/hogehoge/install/build
```

これでインストール作業の準備は完了である．

<a id="GCCのコンパイルとインストール"></a>
<a href="#GCCのコンパイルとインストール"></a>
## GCC のコンパイルとインストール
GCC のコンパイルに必要な依存関係ファイルをダウンロードする．
```
$ cd ../gcc-6.3.0/
$ ./contrib/download_prerequisites
```
makeファイルを作成する  
(ここでは `$HOME/local/gcc-6.3.0` にインストールする設定を与えている)
```
$ cd ../build
$ ../gcc-6.3.0configure --prefix=$HOME/local/gcc-6.3.0 --with-local-prefix=$HOME/local/libgcc63 --enable-checking=release --disable-multilib --enable-languages=c,c++,fortran
```

エラー無くmakefileが生成されたら，コンパイルおよびインストールを行う．
```
$ make
$ make install
```

`$ make` は早くとも30分～1時間程度はかかる．  
無事に `$ make install` が完了したら，`.bashrc` の末尾に以下の設定を追記する．  
あるいは[おまけ2](#複数のgccバージョンの切り替え)のスクリプトを利用する．

```bash
target=$HOME"/local/gcc-6.3.0"  # prefixで指定したGCCのフォルダ
export PATH=${target}/bin:$PATH
export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${target}/lib64:$LD_LIBRARY_PATH
export LD_RUN_PATH=${target}/lib:$LD_RUN_PATH
export LD_RUN_PATH=${target}/lib64:$LD_RUN_PATH
```

.bashrcを変更後，__ログインしなおして__ 新しくインストールしたGCCとパスが有効になっているかどうかを確認する  

```bash
$ gcc -v
Using built-in specs.
COLLECT_GCC=gcc
COLLECT_LTO_WRAPPER=/home/hogehoge/local/gcc-6.3.0/libexec/gcc/x86_64-pc-linux-gnu/6.3.0/lto-wrapper
Target: x86_64-pc-linux-gnu
Configured with: ../gcc-6.3.0/configure --prefix=/home/hogehoge/local/gcc-6.3.0 --with-local-prefix=/home/hogehoge/local/libgcc63 --enable-checking=release --disable-multilib --enable-languages=c,c++,fortran
Thread model: posix
gcc version 6.3.0 (GCC)

$ echo $PATH             # 各々追加したパスが先頭にあるか確認する
$ echo $LD_LIBRARY_PATH
$ echo $LD_RUN_PATH
```
これでGCCのインストールは完了である．
また，`.bashrc` に追記した部分を削除すればデフォルトの環境に戻すこともできる．

<a id="OpenMPIのコンパイルとインストール"></a>
<a href="#OpenMPIのコンパイルとインストール"></a>
## OpenMPI のコンパイルとインストール
既存のbuildディレクトリを使いまわすので gcc 関係のファイルの消去
```
$ rm * -r
```
makeファイルの作成
(ここでは `$HOME/local/openmpi-2.0.1` にインストールする)
```
$ ../openmpi-2.0.1/configure --prefix=$HOME/local/openmpi-2.0.1 --enable-mpi-cxx CC=gcc CXX=g++ F77=gfortran FC=gfortran
```
コンパイルおよびインストール
```
$ make
$ make install
```

GCCと同様に， `.bashrc` に下記のように追記し，パスを登録する．  
あるいは[おまけ2](#複数のgccバージョンの切り替え)のスクリプトを利用する．
```bash
target=$HOME"/local/openmpi-2.0.1"      # prefixで指定したOpenMPIのフォルダ
$ export PATH=${target}/bin:$PATH
$ export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
```
ログインしなおして以下のコマンドの出力を確認する．  
MPIプログラムのコンパイルコマンドが有効かどうか確認
```
$ mpicxx -v
```
出力例： `CXX=[c++ compiler]` で指定したコンパイラのバージョン情報が表示される．  

OpenMPIライブラリのバージョンを確認する
```
$ ompi_info
```
たくさん情報が出てくるが，2行目の Open MPI: (version) や，
後ろのほうの各対応APIのバージョンの右の ", Component v (version)" の数字が
インストールしたものとあっているか確認する．

<a id="FFTwのコンパイルとインストール"></a>
<a href="#FFTwのコンパイルとインストール"></a>
## FFTw のコンパイルとインストール
既存のbuildディレクトリを使いまわすのでOpenMPI関係のファイルの消去
```
$ rm * -r
```

makeファイルの作成  
参考：http://www.fftw.org/doc/Installation-on-Unix.html  
(ここでは `$HOME/local/fftw-3.3.6` にインストールする)  
```
$ ../fftw-3.3.6/configure --prefix=$HOME/local/fftw-3.3.6  \
--enable-mpi --enable-threads --enable-openmp  \
--enable-static --enable-shared  \
--enable-sse2 --enable-avx --with-g77-wrappers  \
--enable-float --enable-sse
```
コンパイルおよびインストール
```
$ make
$ make install
```
このFFTw (32bit単精度版) を使用する場合，コンパイル時に
```bash
-I$HOME/local/fftw-3.3.6/include
-L$HOME/local/fftw-3.3.6/lib
-lfftw3f
-lfftw3f_mpi
-lm
```
をコンパイルオプションに追加する．

FDPSのParticleMesh機能拡張は単精度版のfftwを使用するが，ほかの精度のバージョンも使いたければ
configure 設定の最後につけた `--enable-float --enable-sse` を以下のように書き換える．
```
(消去，このオプションを付けない)  : 標準の倍精度(64bit)浮動小数点
--enable-long-double           : 標準の４倍精度(128bit)浮動小数点
--enable-quad-precision        : 非標準 __float128 四倍精度(128bit)浮動小数点
```
これらの詳細は参考HPを参照のこと．
また，異なる精度のFFTwの指定は上記コンパイルオプションの `-lfftw3*` を書き換えることで行う．上記は32bit単精度版の場合で，64bit倍精度版なら `-lfftw3 -lfftw3_mpi` となる．


GCC, OpenMPIと同様に `.bashrc` に下記のように追記し，パスを登録する．  
あるいは[おまけ2](#複数のgccバージョンの切り替え)のスクリプトを利用する．
```bash
target=$HOME"/local/fftw-3.3.6"      # prefixで指定したFFTwのフォルダ名
$ export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
```
FFTwを使用するプログラムをコンパイルしてできた実行ファイル(例えば a.out )にlddコマンドを使ってFFTwが存在しているかどうかを確認できる．
```
$ ldd a.out
```


<a id="FDPSのparticle_mesh拡張機能のコンパイル"></a>
<a href="#FDPSのparticle_mesh拡張機能のコンパイル"></a>
## おまけ1：FDPSのparticle_mesh拡張機能のコンパイル
ダウンロードしたFDPSのファイルを適当なフォルダ(ここでは `$HOME/local/FDPS-3.0` )に展開しておく．

Makefileの編集(configureは付属していない)  
上記のように環境を構築しているのなら,
```bash
$ cd $HOME/local/FDPS-3.0/src/particle_mesh
$ emacs Makefile
```
このMakefileのコンパイルオプションとMPIライブラリの参照先を編集する．

まずデバッグ用の派生版を作る
```bash
CC = mpicxx
CFLAGS = -O0 -DMPICH_IGNORE_CXX_SEEK -g3 -Wall
INCLUDE_FFTW = -I$HOME/local/fftw-3.3.6/include
```
デバッグ版のコンパイル
```
$ make
```
`libpm.a` (デバッグ版)が作成される．これを別名で (例えば `libpm_debug.a` 等) 保存する．
```
$ mv libpm.a libpm_debug.a
```
Makefile を正式版仕様に再編集
```
$ emacs Makefile
CFLAGS = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK
```
デバッグ用バージョンの中間生成物を除去
```
$ rm -r *.o
```
コンパイル
```
$ make
```
`libpm.a` (正式版)が生成される．


<a id="複数のgccバージョンの切り替え"></a>
<a href="#複数のgccバージョンの切り替え"></a>
## おまけ2: 複数のgccバージョンの切り替え
複数の gcc のバージョンを切り替えたい場合，下記のようなスクリプトを作成しておき 自分の `.bashrc` に読み込むと便利である．

```bash
#!/bin/sh
# additional environment setting for FDPS

#--- select GCC & library (write directory name)
#use_gcc="default"
use_gcc="gcc-6.3.0"

use_ompi="openmpi-2.0.1"
use_fftw="fftw-3.3.6"

#--- user install path
usr_env=$HOME"/local"

#------ symbolic link for default environment
default_env=$HOME"/local_path"

#------------------------------------------------------------
#   All settings in above
#------------------------------------------------------------
#------ PATH for GCC (user installed)
target=${usr_env}'/'${use_gcc}
if [ ${use_gcc} = "default" ]; then
  :
else
  #--- add symbolic link
  export PATH=${default_env}:$PATH

  #--- add user install GCC
  if [ -f ${target}/bin/g++ ]; then
    export PATH=${target}/bin:$PATH
    export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=${target}/lib64:$LD_LIBRARY_PATH
    export LD_RUN_PATH=${target}/lib:$LD_RUN_PATH
    export LD_RUN_PATH=${target}/lib64:$LD_RUN_PATH
  fi
fi

#--- PATH for OpenMPI
target=${usr_env}'/'${use_ompi}
if [ -f ${target}/bin/mpicxx ]; then
  export PATH=${target}/bin:$PATH
  export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
fi

#--- PATH for FFTw
target=${usr_env}'/'${use_fftw}
if [ -f ${target}/include/fftw3.h ]; then
  export LD_LIBRARY_PATH=${target}/lib:$LD_LIBRARY_PATH
fi
```

上記のスクリプトファイルの名称が `user_opt.sh` という名前で自分のホームディレクトリに保存してあるとすると，下記のように `.bashrc` に追記することで設定スクリプトを読み込める
```bash
source $HOME/user_opt.sh
```

一番上の部分の `use_gcc` , `use_ompi` , `use_fftw` に自分がインストールした( `--prefix` に指定した)フォルダ名を，
`usr_env` にそれらのインストールフォルダが存在するディレクトリ (上記の例では `$HOME/local` )を入力するだけでパスの設定が完了する．  
また，違うバージョンを後から追加した場合でも，同じインストールディレクトリにインストールしていれば上記変数のフォルダ名のバージョンを新しくインストールしたものに変えれば適用できる．  
gcc を環境デフォルトのものに戻したくなった場合には，スクリプトの `use_gcc="default"` のコメントアウトを戻し，バージョン指定の方をコメントアウトすればよい．
