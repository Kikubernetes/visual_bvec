# bvec視覚化

本ツールは拡散 MRI における勾配方向（bvec）の空間分布を、直感的かつ定量的に評価する MATLAB ベースのツールである。

bvec の均等性は、拡散指標やトラクトグラフィの精度に影響する重要な要素だが、ベクトル座標の数値だけでは方向分布の特徴を把握しにくい。本ツールは、bvecの軸方向分布を球面上に可視化し、複数の補完的指標を用いて拡散 MRI プロトコルの方向分布の均等性を評価する。

## 特徴

- 可視化するシェルを選択可能（一部のシェルだけ、もしくは全シェルなど）
- 各シェルの軸方向分布を球面上に密度ヒートマップとして可視化
- シェルごとのさまざまな均等性指標による評価をまとめてHTML 形式でレポート出力
- ヒートマップやグラフはPNG、数値はエクセルファイルとして自動保存

### 密度ヒートマップ

球面上で方向密度を可視化する。 
軸数が異なるデータ間でもスケールを揃えることで、偏りの違いを直感的に比較できる。

### 均等性指標

#### 1. Spherical Voronoi

球面上の各方向が担当する領域面積を求め、面積のばらつきを評価する。 
均等な方向分布では、各セル面積は近い値になり、変動計数やmin/maxは小さくなる。

#### 2. Spherical harmonicsによる指標

方向分布を球面調和関数で展開し、低次・高次成分の偏りを定量化する。
分布の全体的な不均一性を捉えるのに役立つ。

#### 3. 最近傍角度

各方向について最も近い方向との角度を計算する。
局所的な密集や疎な領域を把握できる。

### 出力

本ツールでは、以下のような出力が生成される。

- 球面密度ヒートマップ
- Voronoi cell area のヒストグラム
- 最近傍角度のヒストグラム
- 指標一覧テーブル
- 上記をまとめたHTML レポート

### 想定される用途

- 拡散 MRI プロトコル設計時の品質確認
- 既存プロトコル間の方向分布比較
- bvec の subset 作成時の均等性確認

### 動作環境

- MATLAB

### インストール

#### 1 リポジトリを clone

```bash
git clone https://github.com/Kikubernetes/visual_bvec.git
```

#### 2 MATLAB を起動し、このリポジトリをパスに追加 

```matlab
addpath(genpath('path/to/visual_bvec'))
```

## 実行例

まず評価したいbvec / bval ファイルを用意する。FSL形式推奨 （例：dcm2niixの出力ファイルなど）。

### ヒートマップを表示する

bvec, bvalファイルのあるディレクトリにGUIで移動（または絶対パスでファイルを指定してもよい）

以下のように表示に必要な値を入力

```matlab
bvec_file       = 'my_protocol.bvec';
bval_file       = 'my_protocol.bval';

target_b_list   = [1000 2000]; % ヒートマップを表示したいb値をリスト形式で指定
tol             = 50;　　　　　　% b値のゆれを許容する範囲（995、2015等） 通常50でOK
sigma_deg       = 20;          % 大きいほどスムーズな表示になる。通常は20でOK
use_antipodal   = true;　　　　　% 反対方向の軸も使用。通常はtrueを推奨
```

```matlab
% 実行
ss_plot_rel_heatmap( ...
    bvec_file, bval_file, target_b_list, tol, sigma_deg, use_antipodal);
```

ヒートマップは自動的に表示（インタラクティブに回転可能）  
スクリーンショットは以下のような名前でPNGとしてカレントディレクトリに出力される  
`my_protocol_b[1000_2000].png`

### レポートを出力する

表示に必要な値を入力

```matlab
dataset_name  = 'my_protocol';  % レポートに表示される名前。自由に指定できる
bvec_file     = 'my_protocol.bvec';
bval_file     = 'my_protocol.bval';

target_b_list = [1000 2000]; % ヒートマップをみたいb値をリスト形式で指定
tol           = 50;　　　　　　% b値のゆれを許容する範囲（995、2015等） 通常50でOK
lmax          = 8;           % 球面調和関数に使用する次数 
html_name     = 'protocol1'; % フォルダ名になる。複数回実行する場合は上書き注意
```

```matlab
% 実行
bvec_uniformity_report( ...
    dataset_name, bvec_file, bval_file, target_b_list, tol, Lmax, html_name);
```

カレントディレクトリに`html_name`で指定した名前のフォルダができる  
（上記の場合は `protocol1`というフォルダ）  
フォルダ内の `protocol1.html`をダブルクリックするとレポートがブラウザに表示される

#### プロトコルが複数のブロックに分かれている場合

複数ファイルに分かれているbvecを結合して評価したい場合、予め結合済みのbvec, bvalファイルを作成しておく。レポジトリ内に結合用のスクリプトが用意されている。  
DWI1, DWI2の2つに分かれているとすると、

```bash
ombine_bvecs_and_bvals.sh DWI1.bvec DWI2.bvec
```

bvecのみ指定すれば、bvalも自動的に結合される。結合後のファイル名は`combined_DWI1_DWI2.bvec`, `combined_DWI1_DWI2.bval`となる

