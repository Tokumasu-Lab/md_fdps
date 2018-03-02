### note for ./generic_ext

とくに抽象性の高いライブラリについて  
詳細は Doxygen コメントまたは `./unit_test/` にあるテストコードを参照．

#### 利用方法
このフォルダの `molecular_dynamics_ext.hpp` をインクルードする．

#### VEC_EXT
FDPSに付属する PS::Vector3<> の拡張．  
`dot(v1, v2)` : 内積の別名  
`cross(v1, v2)` : 外積の別名  
`rot_x(v, f)`, `rot_y(v, f)`, `rot_z(v, f)` : X, Y, Z各軸周りの回転  
を追加する．

#### MD_EXT::fixed_vector<T,S>
std::vector<T> と同様のインターフェイスで使えるようにラップした std::array<T,S>.  
定義時に指定したサイズまでしか伸長できないがFullParticleなどのMPI通信でやり取りされるデータクラス内にそのまま定義できる．  

#### COMM_TOOL
std::vector<> をはじめとする，いくつかのSTLコンテナとその組み合わせについて，自動的に serialize, communicate, deserialize を行う．  
基本的な集団通信 `broadcast()` , `gather()` , `scatter()` , `allGather()` , `allToAll()` , が使用可能．  
具体的に利用可能なコンテナ，挙動はそれぞれ  
```
./unit_test/gtest_comm_tool_broadcast.cpp
./unit_test/gtest_comm_tool_gather.cpp
./unit_test/gtest_comm_tool_scatter.cpp
./unit_test/gtest_comm_tool_allGather.cpp
./unit_test/gtest_comm_tool_allToAll.cpp
```
を参照．  
より複雑なデータ構造，あるいは通信性能の最適化が必要な場合は [boost.MPI](https://boostjp.github.io/tips/mpi.html) および [boost.Serialization](https://boostjp.github.io/tips/serialize.html) の利用を推奨．

#### STR_TOOL
`split(str, delim)` : `std::string` を指定の区切り文字で分割した `std::vector<std::string>` に変換する．  
`removeCR(str)` : Windowsで編集したテキストファイルの行末に追加される CR コードを除去する．

#### FS_TOOL
ファイル出力用の汎用ツール．  
`make_directory(dir_name, rank)` : ディレクトリ `dir_name` が存在しない場合に新規作成する．  
`FilePrinter` : 出力ファイルの管理用クラス．担当する `rank` からの入力のみをファイルに出力する．  
`file_load(file_name, data_list, rank)` : `file_name` の各行ごとに `std::vector<T>` に読み込んだ `data_list` を作成する．行の読み込みは `void read_line(str)` 関数をデータ型に用意しておく．あるいは `std::vector<std::string>` であれば各行ごとにそのまま格納される．

#### hash_tuple
`std::tuple<>` を `std::unordered_set<>` , `std::unoudered_map<>` のキーとして使用するためのハッシュ関数を提供する．  
`tuple` の各要素のハッシュの生成には `std::hash<>` を用いる．

#### chrono_str
`to_str_h_m_s_ms(chrono)` : `std::chrono` の値を [時間]:[分]:[秒].[ミリ秒] 表記の `std::string` に変換する．  
同様のフォーマットで  
`to_str_m_s_ms(chrono)`  
`to_str_s_ms(chrono)`  
`to_str_ms(chrono)`  
が定義されている．

#### MD_EXT::boltzmann_dist
無次元化マクスウェル・ボルツマン分布生成器を生成する．  
標準の分布生成器と同様に，`[0.0~1.0)` の範囲の実数を引数に与えると無次元化マクスウェル・ボルツマン分布に従う実数を返す．  
内部テーブルの範囲，分解能をコンストラクタ引数または `init()` 関数でカスタマイズ可能．

#### MD_EXT::CellIndex
セルインデックス法によるネイバーリスト生成クラス．MPI通信を行わない．  
系の初期化時の分子配置の衝突判定など，シングルプロセスで実行したい処理でネイバーリストが欲しい場合に用いる．  

#### Normalize
系のサイズ，および正規化空間と実空間の変換関数を提供する．  
FDPSのParticleMesh拡張を用いる場合，FDPSに与える粒子座標は `[0.0~1.0)` である必要があるが，その他実空間の方が都合がいい処理(分子内力の計算など)を行う前後にこれを用いて位置を変換する．

#### MD_EXT::basic_connect
`MD_EXT::fixed_vector<T,S>` のインターフェイスを原子間結合の取り扱い向けに変更したもの．

#### IntraPair
原子間結合の情報から，分子内マスク，angleのペア，dihedralおよびimproper torsionのペアを生成する．  
粒子の持つどの情報を参照するか，という部分を関数オブジェクトでカスタマイズ可能．  
詳細は `./unit_test/gtest_intra_pair.cpp` を参照．
