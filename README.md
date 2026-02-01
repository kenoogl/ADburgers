# ADburgers.jl

ADburgers.jl は、1次元バーガース方程式 (Viscous Burgers' equation) を解くための Julia パッケージです。
高精度な時間積分法（Teylor展開法）と自動微分的なアプローチ（係数の再帰的計算）を用いて、粘性係数や初期条件に応じたシミュレーションを行います。

## 特徴
*   **Taylor Series Method**: 時間方向の高次精度積分（係数の自動計算）
*   **5つのベンチマーク問題**:
    *   **Problem 1**: Wood の厳密解（正弦波の減衰）
    *   **Problem 2**: Cole-Hopf 変換による解（正弦波初期条件）
    *   **Problem 3**: Cole-Hopf 変換による解（放物線初期条件）
    *   **Problem 4**: 衝撃波的な解挙動
    *   **Problem 5**: 区分定数初期条件（不連続性への対応）
*   **検証済み**: 各問題に対する自動化された検証スクリプト完備

## インストール
リポジトリをクローンして、Julia のパッケージモードで `instantiate` してください。

```bash
git clone https://github.com/kenoogl/ADburgers.git
cd ADburgers
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 使用方法

### 基本的なソルバーの実行
```julia
using ADburgers
using Plots

# 問題の取得 (例: Problem 1)
prob = get_problem(1)

# ソルバー設定
# N: 格子点数, dt: 時間刻み, K: Taylor展開次数
spec = SolverSpec(N=40, Δt=0.0001, K=6, save_times=[0.0, 0.5, 1.0])

# 計算実行
sol = solve(prob, spec)

# 結果のプロット (RecipesBase による plot 拡張)
plot(sol, title="Burgers Equation (P1)")
```

### ベンチマーク検証の実行
検証用スクリプトを実行することで、全問題の整合性を確認できます。
```bash
julia --project=. scripts/verify_benchmarks.jl
```
このスクリプトは `plots/` ディレクトリに結果のグラフを生成し、誤差基準（L2/Linfノルムなど）に基づいて合否判定を行います。

## ディレクトリ構成
*   `src/`: ソースコード (`ADburgers.jl`, `core.jl`, `problems.jl` 等)
*   `scripts/`: 検証・比較用スクリプト
*   `test/`: ユニットテスト
*   `.kiro/`: プロジェクト仕様書・要件定義など

## ライセンス
MIT License
