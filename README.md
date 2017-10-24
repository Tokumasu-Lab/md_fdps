# md_fdps


大規模N体シミュレーションフレームワーク[FDPS](https://github.com/FDPS/FDPS)を利用した汎用分子動力学シミュレーションコード開発ベース……の予定

### 動作環境
**C++11規格を使用**

開発に使用している環境は
  - GCC 6.3.0
  - OpenMPI 2.1.1
  - FFTw 3.3.6  

である．

開発環境の構築手順については `./ENV-NOTE.md` を参照  
`makefile` のインクルードパス，ライブラリのリンク指定は上記で構築した環境への静的リンクを用いているので適宜変更する．

### コンパイル&実行
原子モデルの識別に用いるヘッダファイル `./src/enum_model.hpp` は  
`./model/` フォルダにある `*.mol2` ファイルを用いて  
```
$ ./script/convert_model_indicator.py
```
で生成する．その後に  
```
$ make
```
で実行ファイル `md_fdps.x` が生成される．

計算条件は  
```
./condition_molecule.imp
./condition_sequence.imp
```
にそれぞれ記入する．

準備が整ったら
```
$ mpirun -n [n_proc] ./md_fdps.x
```
で実行する． [n_proc] は任意のMPIプロセス数で，ParticleMeshの仕様により2以上である必要がある

#### Hybrid並列実行について
コンパイルオプションに `-DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp` を追加する．(MPIは必須)

OpenMPI-2.1.1 + OpenMP の環境で，Xeon E5-2690v3 (12core) を 2CPU のノード4台で各CPUにMPIプロセスを1つ配置し，CPU内部ではOpenMP並列化をおこなう例
```
$ export OMP_NUM_THREADS=12
$ mpirun -n 8 -bind-to socket -npersocket 1 md_fdps.x
```

### 実装目標と現状
  - モデル，設定の読み込み
    - Ar, H2Oは付属 ( `./model/` )．任意に追加可能．モデルを増やした際には `./script/convert_model_indicator.py` を再実行して再コンパイル
  - MPI通信の追加ラッパー
    - broadcast 用のSTLコンテナアダプタ．使い方は `./unit_test/comm_tool.cpp` を参照
  - Intra相互作用ペアのマネージャ **(要改良)**
  - 系の初期化
  - 中断データの生成 **(未実装)** ，続行 **(未実装)**
  - 基本的な古典相互作用
    - LJ, Coulomb, Bond, Angle, dihedral torsion **(仮実装)**, improper torsion **(未検証)**
  - 基本的な時間積分
    - velocity verlet
  - 基本的な拡張系制御
    - NVT, NPT **(検証中)**
  - 基本的な解析
    - RDF **(未実装)**, MSD **(未実装)**
  - 基本的な可視化
    - VMD

コメントを[Doxygen](http://www.doxygen.jp)に対応するフォーマットに変更中

### VMDによる動画作成
`md_fdps.x` を実行し，`./pdb` ， `./posdata` が生成されたものとする．  
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


### Contact
東北大学　流体科学研究所  
徳増研究室 md_fdps 開発チーム  
contact.md-fdps -@- nanoint.ifs.tohoku.ac.jp  
( " -@- " を "@" に置き換えてください)
