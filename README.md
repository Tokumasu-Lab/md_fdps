# md_fdps

大規模N体シミュレーションフレームワーク[FDPS](https://github.com/FDPS/FDPS)を利用した汎用の分子動力学シミュレータ実装  

### 概要
基本的な分子動力学解析の流れ：  

 - 分子モデルの読み込み  
 - 初期構造の作成  
 - シミュレーションの実行  
 - 出力結果の解析

までの一通りの処理を実行可能なシンプルなシステムを提供する．

そのままでも利用は可能だが，ユーザがオリジナルの相互作用や解析処理を追加することで，  
少ない手間で高度に並列化された新規の分子動力学シミュレータを構築する基盤として用いることを想定している．

### 動作環境
**C++11規格を使用**

開発に使用している環境は
  - FDPS 4.0c
  - GCC 6.3.0
  - OpenMPI 2.1.2
  - FFTw 3.3.7  

である．

開発環境の構築手順については `./ENV-NOTE.md` を参照  
`makefile` のインクルードパス，ライブラリのリンク指定は上記で構築した環境への静的リンクを用いているので適宜変更する．  
あるいは下記のように環境変数を追加する．
```
export FDPS_ROOT=[FDPSの src フォルダが見えるディレクトリ]
export FFTW_ROOT=[FFTwの include フォルダが見えるディレクトリ]
```

### コンパイル&実行
開発環境の構築ができたら，下記のスクリプトに実行権限を付与して `make` を実行する．
```
$ chmod ugo+x ./script/convert_model_indicator.py
$ make
```
初期化処理用の実行ファイル `md_init.x` と  
時間積分用の実行ファイル `md_fdps.x` が生成される．

計算条件は  
```
./condition_molecule.inp
./condition_sequence.inp
```
にそれぞれ記入する．  
準備が整ったら
```
$ ./md_init.x
```
で系の初期状態を生成した後，
```
$ mpirun -np [n_proc] ./md_fdps.x
```
で分子動力学シミュレーションを実行する．  
ここで， `[n_proc]` は任意のMPIプロセス数で，FDPSのParticleMesh拡張の仕様により2以上である必要がある

#### Hybrid並列実行について
コンパイルオプションに `-DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp` を追加する．  
(先述の通り，MPIは必須)

OpenMPI-2.1.2 + OpenMP のHybrid並列化で，Xeon E5-2690v3 (12core) を 2CPU のノード4台 (合計8CPU) で各CPUにMPIプロセスを1つ配置し，CPU内部ではOpenMP並列化 (各12スレッド) をおこなう場合の実行コマンドは下記のようになる．
```
$ mpirun -np 8 --bind-to socket -npersocket 1 -x OMP_NUM_THREADS=12 md_fdps.x
```
現状の実装では，メモリ容量や通信量の削減の必要が出ない限りHybrid化せず flat-MPI で実行したほうが速い．

### 実装目標と現状
  - モデル，設定の読み込み
    - All-Atom, flexibleモデルのAr, 水，プロパノール2種，トルエンは付属 ( `./model/` )．
    - 分子モデルは任意に追加可能．
    - 分子モデルを追加した際には `./script/convert_model_indicator.py` を実行して  
    `./src/enum_model.hpp` を更新してからコードをコンパイルする．
  - Intra相互作用ペアのマネージャ
    - 分子間相互作用に対する分子内 mask の係数および範囲はモデルごとに任意に設定可能．
    - 具体的な挙動は Doxygenコメントおよび `./unit_test/gtest_intra_pair.cpp` を参照．
  - 系の初期状態の作成
    - 初期状態の作成は `md_init.x` で行い，出力された resume データを `md_fdps.x` で読み込む．
    - `./condition_molecule.inp` で指定された分子を指定の個数，指定の分子間距離でランダムに配置し，  
     `./condition_sequence.inp` の 0 step 目に相当する温度のマクスウェル・ボルツマン速度分布を与える．
  - resume データの出力およびシミュレーションの再開
    - ascii形式単一ファイル入出力, 浮動小数点は丸め処理を避けるため16進数表記
  - 基本的な古典相互作用
    - 分子間：
      - LJ (6,12形)
      - Coulomb
    - 分子内：
      - Bond
      - Angle
      - dihedral torsion
      - improper torsion
    - 分子間力の特定のペアに任意の係数の分子内マスクを考慮
  - 基本的な時間積分
    - velocity verlet
  - 基本的な拡張系制御
    - NVT
    - NPT
  - 基本的な解析
    - Radial Distribution Function (RDF) **(ユニットテストまで実装)**
    - Mean Square Displacement (MSD) **(ユニットテストまで実装)**
    - Clustering **(ユニットテストまで実装)**
  - 基本的な可視化
    - VMD

一部のコメントで[Doxygen](http://www.doxygen.jp)に対応

### VMDによる動画作成
`md_fdps.x` を実行し，`./pdb` と `./posdata` が生成されたものとする．  
ここで，  
```
$ ./script/VMDmovie_convert.py
```
を実行し，生成された  
```
./vmd_movie.pdb
./vmd_movie.crd
```
を[VMD](http://www.ks.uiuc.edu/Research/vmd/)で読み込み可視化する．


### Unit test について
`./unit_Test/gtest_***.cpp` の名前のテストコードでは [Google C++ Test framework](https://github.com/google/googletest) を使用している．  
現状は version 1.8.0 を使用．

付属の `makefile` あるいは実行スクリプト `./gtest_run.sh` を用いてユニットテストをコンパイル・実行する場合，環境変数として以下を設定する必要がある．
```bash
export GTEST_ROOT=[Google C++ Test framework 本体の解凍後に cmake, make を実行したディレクトリ]
```
また，テストスクリプト実行時に以下のように設定を変更できる．
```bash
$ ./gtest_run.sh -j [コンパイル時の並列数] -m [MPI並列数] -o [OpenMP並列数]
```

そのほか，分子モデルや計算条件の設定ファイルに関する確認は以下のテストで行う．  
こちらは google test を使用していないので，`md_fdps.x` 本体のコンパイルに必要な環境と同じ環境のみでテスト可能である．
```bash
#--- 計算条件の読み込み
$ make test_condition
#--- 分子モデルパラメータの読み込みテスト
$ make test_model
#--- 分子内力パラメータの定義エラーテスト
$ make test_param
```

各ユニットテストの実行ファイルは `./test_bin/` 以下に生成される．
上記のテストのうち，分子モデルのテストは引数に任意の分子モデル名を与えることで，あとから追加した分子モデルのテストも可能である．
```bash
./test_bin/test_model [model_name]

./test_bin/test_param [model_name]
```
ここで， `[model_name]` は `***.mol2` , `***.param` ファイルの拡張子を除いたファイル名の部分で，両方のパラメータファイルがセットになっている必要がある．

### 機能実験用の定義済マクロ

| マクロ | 効果 |
|:-------|:-----|
| REUSE_INTERACTION_LIST | FDPSの相互作用リストと分子内ペアリストを毎ステップ作り直さずに cycle_dinfo 回使いまわす．現状の PS::ParticleMesh の仕様の制限から通信量が増加する．__動作検証中__  |
| FORCE_NAIVE_IMPL | 分子内マスクを短距離相互作用カーネル内で直接評価する．桁落ちの心配はなくなるが性能が低下する．(デバッグ用) |
| CHECK_FORCE_STRENGTH | 時間積分が明らかに破綻するような巨大な力が働いた粒子を検知しレポートする．(デバッグ用) |
| DEBUG_COMM_TOOL | std::vector や std::pair の入れ子など，複雑な形状のデータを COMM_TOOL の集団通信に渡した際に再帰呼び出しが正しいか確認するための情報を出力する．(デバッグ用) |


### Contact
東北大学　流体科学研究所  
徳増研究室 md_fdps 開発チーム  
contact.md-fdps -@- nanoint.ifs.tohoku.ac.jp  
( " -@- " を "@" に置き換えてください)
