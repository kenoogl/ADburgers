# 基本設計: adburgers

## 1. モジュール構成（単一モジュール + ファイル分割）

```text
src/
├── ADburgers.jl
├── types.jl            # ProblemSpec, SolverSpec, Solution, ReferenceSpec, TaylorArrays
├── core.jl             # solve, step!, taylor_coeff!, horner_update!
├── problems.jl         # Problem 1-5 (u0, bc, exact/ref)
├── factory.jl          # get_problem(id)::(ProblemSpec, ReferenceSpec, SolverHint)
├── analysis.jl         # errors, ref_self_convergence, check_acceptance(ACC-*)
└── visualization.jl
scripts/
└── verify_benchmarks.jl
test/
└── runtests.jl
```

------

## 2. データ構造詳細（types.jl）

### 2.1 型（境界関数のシグネチャ固定）

```julia
const BCTime = Function   # (t::Float64) -> Float64
const BCCoef = Function   # (k::Int, t::Float64) -> Float64
```

### 2.2 ProblemSpec

```julia
struct ProblemSpec
    a::Float64; b::Float64; ν::Float64
    t0::Float64; tmax::Float64
    u0::Function            # (x::Float64) -> Float64
    bc_left::BCTime;  bc_right::BCTime
    bc_left_coeff::BCCoef; bc_right_coeff::BCCoef
end
```

### 2.3 ReferenceSpec（P2/P3 の核心）

```julia
struct ReferenceSpec
    has_exact::Bool           # P1,P4 true
    has_reference::Bool       # P2,P3 true
    Nfourier::Int             # default 1000, min 500
    rtol::Float64             # default 1e-12
    atol::Float64             # default 1e-14
    enforce_self_convergence::Bool  # P2,P3 true
    self_conv_tol::Float64    # 1e-10
end
```

### 2.4 SolverSpec（save_timesは格子一致必須）

```julia
struct SolverSpec
    N::Int; Δt::Float64; K::Int
    save_times::Vector{Float64}  # t0..tmax 内、格子一致必須
end
```

### 2.5 Solution（返り値固定）

```julia
struct Solution
    x::Vector{Float64}
    t::Vector{Float64}         # = save_times（検証後）
    u::Matrix{Float64}         # size (N+1, Nt)
    meta::Dict{String,Any}
end
```

### 2.6 TaylorArrays（U規約：正式決定）

> **本設計書の正式規約：**
> `U[i+1, k+1] == (u_i)_k`、`U` の shape は **(N+1, K+1)**

```julia
struct TaylorArrays
    U::Matrix{Float64}  # (N+1, K+1)
    T1::Vector{Float64} # (N+1)
    T2::Vector{Float64} # (N+1)
    T3::Vector{Float64} # (N+1)
end
```

------

## 3. 前処理・入力検証（core.jl: solve 冒頭）

### 3.1 save_times 格子一致判定（tolを定義）

- `m = (t - t0) / Δt`
- `abs(m - round(m)) ≤ tol_m` を満たす必要
- **設計で固定する tol：**
  - `tol_m = 1000 * eps(Float64) * max(1.0, abs(m))`
- 条件に違反したら `ArgumentError`

### 3.2 save_times 範囲チェック

- `minimum(save_times) == t0` を推奨（含まれない場合は自動で先頭に追加するか、エラーにするかを選ぶ）
  - **本設計は明確性優先でエラー**：`t0` を含めることを必須
- `maximum(save_times) ≤ tmax`（超える場合はエラー）

### 3.3 初期条件と境界の整合性チェック（追加）

- `abs(u0(a) - bc_left(t0)) ≤ tol_bc`
- `abs(u0(b) - bc_right(t0)) ≤ tol_bc`
- `tol_bc = 1e-12`（固定）
- 破る場合：
  - デフォルトは **警告＋境界優先で上書き**（metaに記録）
  - オプションで strict モードならエラー（将来拡張）

------

## 4. アルゴリズム設計（core.jl）

### 4.1 計算順序（Rev.3要件を設計に固定）

1. `U[:,1]`（k=0）に現在状態 `u(t)` を保持
2. `k=0..K-1` で Taylor係数を再帰生成（境界係数上書きを先に行う）
3. Horner法で Taylor 和を計算し、`U[:,1]` を `u(t+Δt)` に更新
4. 境界 `u(a,t+Δt), u(b,t+Δt)` を厳密適用
5. 安定性指標更新（CFL/diffusion）
6. NaN/Inf 検出で停止

### 4.2 境界近傍の扱い（レビュー指摘2の明文化）

- 内点ループは **Julia添字 iJ=2..N**（物理 i=1..N-1）
- 中心差分では `U[iJ-1,*]` と `U[iJ+1,*]` を参照するため、境界 `iJ=1,N+1` の係数が **当該kで事前にセット済み**であることが前提。

### 4.3 T2 畳み込みのループ順序（レビュー指摘1の対策）

- **最内ループは i（空間）**にしてSIMDが効きやすい形を採用（Kが小さくても明文化する）
- 設計上の基本順序：`k -> i -> j`

### 4.4 疑似コード：Taylor係数生成（taylor_coeff!）

（**低アロケーション版**、T1/T3はkで一度、T2は iごとacc）

```text
for k = 0..K-1:
  set boundary coefficients for order k:
    U[1, k+1]   = bc_left_coeff(k, t)
    U[N+1,k+1]  = bc_right_coeff(k, t)

  compute T1_k and T3_k for iJ=2..N:
    T1[iJ] = (U[iJ+1,k+1] - U[iJ-1,k+1]) / (2h)
    T3[iJ] = (U[iJ-1,k+1] - 2U[iJ,k+1] + U[iJ+1,k+1]) / (h^2)

  compute T2_k for iJ=2..N:
    acc = 0
    for j = 0..k:
      T1_{k-j}(iJ) = (U[iJ+1,(k-j)+1] - U[iJ-1,(k-j)+1])/(2h)
      acc += U[iJ, j+1] * T1_{k-j}(iJ)
    T2[iJ] = acc

  update next coefficient (order k+1) for iJ=2..N:
    U[iJ,(k+1)+1] = ( -T2[iJ] + ν*T3[iJ] ) / (k+1)
```

（※ここでは `T1_{k-j}` を保存しない方針。将来Kを上げる場合は `T1coef` 追加で高速化可能）

### 4.5 疑似コード：Horner法による Taylor 和（horner_update!）

```text
for iJ = 2..N:
  s = U[iJ, K+1]
  for kk = K..0:
    s = s*Δt + U[iJ, kk+1]
  U[iJ, 1] = s
set boundary:
  U[1,1]   = bc_left(t+Δt)
  U[N+1,1] = bc_right(t+Δt)
```

### 4.6 NaN/Inf 検出後の挙動（明確化）

- `isfinite` を全内点（あるいは全点）でチェック
- 検出したら：
  - `meta["nan_detected"]=true`
  - その時点で `DomainError` もしくは `ErrorException` を投げて停止
    （設計としては **例外停止を標準**とする）

------

## 5. 安定性指標（meta）定義（レビュー指摘4）

### 5.1 指標の定義（固定）

- `CFL(t) = max_i |u_i(t)| * Δt / h`
- `diffusion = ν * Δt / h^2`

### 5.2 記録タイミング

- 各ステップで算出し、全ステップの最大値を
  - `meta["CFL_max"]`
  - `meta["diffusion_max"]`
    に保存。

------

## 6. Problem定義（自己完結化）

> レビュー「Problem具体値が本文にない」を解消するため、**factory.jl に具体値を内包**し、PDF参照不要にする（ただし論文一致の確認は別途）。

### 6.1 get_problem API

```julia
get_problem(id::Int) -> (prob::ProblemSpec, ref::ReferenceSpec, hint::Dict{String,Any})
```

### 6.2 hint（例）

- `"N"`, `"dt"`, `"K"`, `"tmax"`, `"save_times"` を必ず含める

### 6.3 Problem 4（t0=1）の波及（レビュー指摘3）

- `t0` を一般値として、格子一致判定は常に `m=(t-t0)/dt` を用いる
- `test` で **t0≠0 のケース**を必ず含める（設計に明文化）

------

## 7. Cole–Hopf / Fourier 参照解設計（Problem2/3）

レビュー指摘の「積分区間・正規化定数」を設計に追加。

### 7.1 参照解モジュールの責務

- `reference_u(prob_id, x, t, refspec)` を提供
- 係数計算（数値積分）と級数打切り、自己収束判定を含む

### 7.2 数値積分仕様（区間を明示）

- 積分区間：**問題の定義域 `[a,b]`**（一般形にして自己完結）
- 積分は `QuadGK` 等を想定し、`rtol`, `atol` は `refspec` に従う
- 正規化定数・係数式（具体式）は problems.jl に **式コメント＋参照**として明記（実装者が誤らないよう、式をそのまま残す）

### 7.3 参照解自己収束判定（必須）

- `ref_self_convergence(prob_id, t, xs, refspec)` を実装
- 条件：
  [
  \max_{x\in xs} |u_{N_f}(x,t) - u_{2N_f}(x,t)| \le \text{self_conv_tol}
  ]
- 不合格なら：ACC判定の前に停止（参照解が信頼できないため）

------

## 8. 受入基準（ACC）— P2/P3 を追加

レビュー指摘「ACC-P2/P3欠落」を補完。

### ACC-P1-T1（既存）

- 条件：hint のデフォルト設定で実行
- 判定：指定点の max 絶対誤差 ≤ `5e-7`

### ACC-P2-T*（新規）

- 前提：参照解自己収束が **合格**（`ref_self_convergence == true`）
- 判定：指定時刻（hint）で
  - `L∞(u_num - u_ref) ≤ 5e-6`
  - （必要なら Table点比較も追加）

### ACC-P3-T*（新規）

- 前提：参照解自己収束が **合格**
- 判定：指定時刻（hint）で
  - `L∞ ≤ 5e-6`

### ACC-P4-T8（既存）

- `L∞ ≤ 1e-6`

### ACC-P5-FIG（既存）

- 境界保持 + overshoot制限 + NaN/Inf無し

------

## 9. analysis.jl API（固定）

```julia
errors(sol, prob_id, prob, refspec; t)::NamedTuple
check_acceptance(sol, prob_id, prob, refspec)::Dict{String,Bool}
ref_self_convergence(prob_id, t, xs, refspec)::Bool
```

------

## 10. 追加のテスト要件（レビュー反映）

- `t0≠0`（Problem4）の save_times 格子一致テストを必須化
- `u0` と境界不整合時の動作（警告＋上書き）テスト
- NaN/Inf 検出時に例外で止まることのテスト
