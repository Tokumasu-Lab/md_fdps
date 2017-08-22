#md_fdps
***

大規模N体シミュレーションフレームワーク[FDPS](https://github.com/FDPS/FDPS)を利用した汎用分子動力学シミュレーションコード開発ベース……の予定

###C++11を利用
開発環境の構築手順については `./development_environment.md` を参照

###コンパイル&実行
原子モデルの識別に用いるヘッダファイル `./src/md_fdps_enum_model.hpp` は  
`./model/` フォルダにある `*.mol2` ファイルを用いて  

`$ ./script/convert_model_indicator.py`

で生成する．その後に  
`$ make`  
で実行ファイル `md_fdps.x` が生成される．

計算条件は  
`./condition_molecule.imp`  
`./condition_sequence.imp`  
にそれぞれ記入する．

###実装目標と現状
  - モデル・設定の読み込み
    - Ar, H2Oは付属 ( `./model/` )．任意に追加可能
  - MPI通信の追加ラッパー
    - broadcast 用のSTLコンテナアダプタ．詳細は `./unit_test/comm_tool.cpp` を参照
  - Intra相互作用ペアのマネージャ **(要改良)**
  - 系の初期化
  - 中断データの生成，続行 **(未実装)**
  - 基本的な古典相互作用
    - LJ, Coulomb, Bond, Angle, Torsion **(仮実装)**
  - 基本的な時間積分
    - velocity verlet
  - 基本的な拡張系制御
    - NVT, NPT **(検証中)**
  - 基本的な解析
    - RDF **(未実装)**, MSD **(未実装)**
  - 基本的な可視化
    - VMD


###VMDによる動画作成
`md_fdps.x` を実行し，`./pdb` ， `./posdata` が生成されたものとする．  
ここで，  

`$ ./script/VMDmovie_convert.py`  

を実行し，生成された  
`./vmd_movie.pdb`  
`./vmd_movie.crd`  
を[VMD](http://www.ks.uiuc.edu/Research/vmd/)で読み込み可視化する．
