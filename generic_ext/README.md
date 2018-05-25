### note for ./generic_ext

とくに抽象性の高いライブラリについて  
詳細は Doxygen コメントまたは `./unit_test/` にあるテストコードを参照．

#### 利用方法
このフォルダの `molecular_dynamics_ext.hpp` をインクルードする．

#### VEC_EXT
FDPSに付属する PS::Vector3\<T\> の拡張．  
`dot(v1, v2)` : 内積の別名  
`cross(v1, v2)` : 外積の別名  
`rot_x(v, f)`, `rot_y(v, f)`, `rot_z(v, f)` : X, Y, Z各軸周りの回転(ラジアン)  
を追加する．

#### MD_EXT::fixed_vector\<T, S\>
std::vector\<T\> と同様のインターフェイスで使えるようにラップした std::array\<T,S\>.  
定義時に指定したサイズまでしか伸長できないがFullParticleなどのMPI通信でやり取りされるデータクラス内にそのまま定義できる．
詳細な挙動は `./unit_test/gtest_fixed_vector.cpp` を参照．

[boost.Container](http://www.boost.org/doc/libs/1_66_0/doc/html/container.html) が利用可能なら `boost::container::static_vector` の利用を推奨．

#### MD_EXT::logspace_array\<T\>
対数軸上で等間隔な集計を行うための配列．  
実数のインデックスでアクセスし，その境界は等比数列となっている．
以下のように初期化した場合，`[1.0 ~ 100.0)` の範囲を `10` 分割した領域として管理される．
```c++
MD_EXT::logspace_array<int> ls_array;
ls_array.init(1.0, 100.0, 10);
```

インターフェイスや挙動の詳細は `./unit_test/gtest_logspace_array.cpp` を参照．

__注: NumPyやMATLABの logspace(first, last, size) では first と last を10を底とする指数で指定するが，logspace_arrayでは実数で直接与える.__

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
ファイル入出力用の簡易ツール．  
`make_directory(dir_name, rank)` : ディレクトリ `dir_name` が存在しない場合に新規作成する．  
`FilePrinter` : 出力ファイルの管理用クラス．担当する `rank` からの入力のみをファイルに出力する．  
`file_load(file_name, data_list, rank)` : `file_name` の各行ごとに `std::vector<T>` に読み込んだ `data_list` を作成する．行の読み込みは `void read_line(str)` 関数をデータ型に用意しておく．あるいは `std::vector<std::string>` であれば各行ごとにそのまま格納される．

#### hash_tuple::hash_func\<T\>
`std::tuple<T1,T2,...Tn>` を `std::unordered_set<>` , `std::unoudered_map<>` のキーとして使用するためのハッシュ関数を提供する．  
`tuple` の各要素のハッシュの生成には `std::hash<T>` を用いる．

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
詳細な挙動は `./unit_test/gtest_blz_dist.cpp` を参照．

#### MD_EXT::CellIndex\<T\>
cell index法を用いてネイバーリストを作成する．MPI通信を行わない．  
系の初期化時における分子配置の衝突判定など，シングルプロセスで実行したい処理でネイバーリストが欲しい場合に用いる．  
詳細な挙動は `./unit_test/gtest_cell_index.cpp` を参照．

#### Normalize
系のサイズ，および正規化空間と実空間の変換関数を提供する．  
FDPSのParticleMesh拡張を用いる場合，FDPSに与える粒子座標は `[0.0~1.0)` である必要があるが，その他実空間の方が都合がいい処理(分子内力の計算など)を行う前後にこれを用いて位置を変換する．

#### MD_EXT::basic_connect\<T,S\>
`MD_EXT::fixed_vector<T,S>` のインターフェイスを原子間結合の取り扱い向けに変更したもの．
詳細な挙動は `./unit_test/gtest_basic_connect.cpp` を参照．

#### IntraPair
原子間結合の情報から，分子内マスク，angleのペア，dihedralおよびimproper torsionのペアを生成する．  
粒子の持つどの情報を参照するか，という部分を関数オブジェクトでカスタマイズ可能．  
詳細な挙動は `./unit_test/gtest_intra_pair.cpp` を参照．
