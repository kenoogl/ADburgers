# 要件定義: adburgers

# Burgers方程式

## Taylor級数（自動微分的再帰）ソルバー 要件定義（Rev.3／完成版）

------

## 0. 改訂履歴

### Rev.3（本版）

- 数式表記を **LaTeXで完全整形**（添字・符号の曖昧さを排除）
- Taylor再帰式の **符号規約・更新順序を完全固定**
- **受入基準をテストケースID付きで定義**
- Problem 2/3 の **参照解（Cole–Hopf/Fourier）の精度保証要件を追加**
- `solve` API、返り値構造、`save_times` の挙動を一意化
- Taylor係数配列 `U` の **数学添字 ↔ Julia実装添字の対応を明文化**

------

## 1. 目的・スコープ

### 1.1 目的

論文
**“Numerical solution of the Burgers' equation by automatic differentiation”**
で提案された
**Method of Lines + Taylor級数（自動微分的再帰）法**を Julia で忠実に実装し、
Problem 1–5 の **数値表・図を追試再現**する。

### 1.2 対象方程式

1 次元粘性 Burgers 方程式：
$$
\frac{\partial u}{\partial t}

u \frac{\partial u}{\partial x}
= \nu \frac{\partial^2 u}{\partial x^2},
\qquad a < x < b,; t > 0
$$
境界条件：Dirichlet
（Problem 5 のみ非同次）

------

## 2. 入出力仕様

### 2.1 入力：ProblemSpec

| 項目                  | 型         | 説明                           |
| --------------------- | ---------- | ------------------------------ |
| `a, b`                | `Float64`  | 空間領域                       |
| `ν`                   | `Float64`  | 粘性係数                       |
| `t0`                  | `Float64`  | 初期時刻（Problem4 は `t0=1`） |
| `tmax`                | `Float64`  | 最終時刻                       |
| `u0(x)`               | `Function` | 初期条件                       |
| `bc_left(t)`          | `Function` | 左境界                         |
| `bc_right(t)`         | `Function` | 右境界                         |
| `bc_left_coeff(k,t)`  | `Function` | 左境界の Taylor 係数           |
| `bc_right_coeff(k,t)` | `Function` | 右境界の Taylor 係数           |

#### 境界条件規約（必須）

- **内部では必ず `bc_\*_coeff(k,t)` を使用**

- 定数 Dirichlet の場合
  $$
  u(a,t)=C \Rightarrow
  \begin{cases}
  (u)_0 = C \
  (u)_k = 0 \quad (k \ge 1)
  \end{cases}
  $$
  

------

### 2.2 入力：SolverSpec

| 項目         | 型                | 説明                     |
| ------------ | ----------------- | ------------------------ |
| `N`          | `Int`             | 分割数（格子点数 = N+1） |
| `Δt`         | `Float64`         | 時間刻み                 |
| `K`          | `Int`             | Taylor 次数（推奨 1–6）  |
| `save_times` | `Vector{Float64}` | 保存時刻（※制約あり）    |
| `dtype`      | `DataType`        | `Float64` 固定           |

#### `save_times` 規約

- **時間格子と一致する値のみ許可**
- 一致しない場合はエラーを返す（補間は禁止）

------

### 2.3 出力：Solution

```julia
struct Solution
    x    :: Vector{Float64}     # 空間格子
    t    :: Vector{Float64}     # 保存時刻
    u    :: Matrix{Float64}     # size = (N+1, Nt)
    meta :: Dict{String,Any}    # 条件・警告・安定性指標
end
```

------

## 3. 数値手法仕様（コア）

### 3.1 空間離散（Method of Lines）

内部点 ( i = 1,\dots,N-1 )：
$$
\frac{d u_i}{dt}

u_i \frac{u_{i+1}-u_{i-1}}{2h}

\nu \frac{u_{i-1}-2u_i+u_{i+1}}{h^2},
  \quad h = \frac{b-a}{N}
$$


境界点 ( i=0,N ) は Dirichlet により拘束。

------

### 3.2 Taylor 係数の定義

$$
(u_i)_k

\frac{1}{k!}
\left.\frac{d^k u_i}{dt^k}\right|_{t=t_n}
$$



時間更新：
$$
u_i(t_{n+1})

\sum_{k=0}^{K} (u_i)_k (\Delta t)^k
$$


------

### 3.3 再帰式（符号規約：固定）

#### 補助量

$$
(T1_i)_k
\frac{(u_{i+1})*k - (u*{i-1})_k}{2h}
\\
(T2_i)_k

\sum_{j=0}^{k}
(u_i)*j (T1_i)*{k-j}
\\
(T3_i)_k
\frac{(u_{i-1})_k - 2(u_i)*k + (u*{i+1})_k}{h^2}
$$



#### 更新式（**最重要仕様**）

$$
(u_i)_{k+1}
\frac{ - (T2_i)_k + \nu (T3_i)_k }{k+1}
$$



※ 移流項は **必ずマイナス**。

------

### 3.4 計算順序（必須）

1. 各ステップ開始時：`U[:, 1]` (k=0) を設定
2. **k = 0 → K-1** の順で
   - 境界係数を上書き
   - `(T1)_k → (T2)_k → (T3)_k`
   - `(u)_{k+1}` を計算
3. Taylor 和で `u_next` を更新
4. 境界値を再適用

------

## 4. データ構造規約（重要）

### Taylor 係数配列

- 数学添字：((u_i)_k)
- Julia 実装：

```julia
U[i+1, k+1] == (u_i)_k
```

- `U :: Matrix{Float64}` size = `(N+1, K+1)`
- 境界点 `i=0,N` は **各 k で必ず上書き**

------

## 5. 機能要件

### F1. Problem Registry

- Problem 1–5 のプリセット定義
- 推奨 `(N, Δt, K, tmax)` を明示
- Problem 4 は `t0=1` 対応

### F2. Taylor Engine

- 再帰係数計算
- O(K²) 畳み込み（K≲6 前提）

### F3. Time Stepper

- Taylor 和による更新
- 境界再適用

### F4. Exact Solution

- P1, P4：閉形式解
- P2, P3：Cole–Hopf Fourier 級数
- P5：参照解なし（物理挙動評価）

### F5. Analytics

- 表出力
- L∞, L2 誤差
- 収束テスト（Δt/2）

### F6. Visualization

- u(x,t) プロット
- Problem 5 は遷移・ショック形状重視

------

## 6. 参照解（Problem 2/3）精度保証

### Fourier 級数要件

- `Nfourier ≥ 500`（推奨 1000）
- 数値積分：`rtol ≤ 1e-12`, `atol ≤ 1e-14`

### 自己収束チェック（必須）

$$
\max_x |u_{N_f}(x,t) - u_{2N_f}(x,t)| \le 10^{-10}
$$



満たさない場合は **参照解不合格**。

------

## 7. 受入基準（Acceptance Criteria）

### ACC-P1-T1

- Problem1

- 条件：論文 Table 対応条件

- 判定：指定 x 点で
  $$
  |u_{\text{num}} - u_{\text{exact}}| \le 5 \times 10^{-7}
  $$
  

### ACC-P4-T8

- Problem4

- 判定：
  $$
  |u_{\text{num}} - u_{\text{exact}}|_{\infty} \le 10^{-6}
  $$
  

### ACC-P5-FIG

- Problem5
- 判定条件：
  - 境界保持：`u(0)=1`, `u(b)=0`
  - `max(u) ≤ 1 + 10^{-3}`
  - 非物理振動（NaN/Inf）なし

------

## 8. 非機能要件

### 安定性ログ

- `max|u| Δt / h`
- `ν Δt / h²`
- NaN/Inf 検出で即停止

### 性能

- N≤200, K≤6 で実用時間

### 再現性

- 条件・参照解設定を JSON/TOML 保存

------

## 9. テスト要件

### Unit Test

- 境界条件保持
- k=1 が MOL RHS と一致

### Regression Test

- ACC-P1-T1
- ACC-P4-T8
- ACC-P5-FIG

------

## 10. ベンチマーク問題定義 (Problem 1–5)

本プロジェクトでは以下の5つの標準問題を検証対象とする。

**共通支配方程式**:
$$
u_t + u u_x = \nu u_{xx}
$$

### Problem 1: Exact Solution (Wood)

- **領域**: $0 \le x \le 1$
- **パラメータ**: $a=2$ (または $a>1$)
- **初期条件**:
  $$
  u(x,0) = \frac{2\nu\pi\sin(\pi x)}{a+\cos(\pi x)}
  $$
- **境界条件**: $u(0,t)=u(1,t)=0$ (Homogeneous Dirichlet)
- **厳密解**:
  $$
  u(x,t) = \frac{2\nu\pi e^{-\pi^2 \nu t}\sin(\pi x)}{a+e^{-\pi^2 \nu t}\cos(\pi x)}
  $$

### Problem 2: Cole-Hopf / Fourier Series

- **領域**: $0 \le x \le 1$
- **初期条件**:
  $$
  u(x,0) = \sin(\pi x)
  $$
- **境界条件**: $u(0,t)=u(1,t)=0$
- **厳密解**:
  $$
  u(x,t) = \frac{2\pi\nu \sum_{n=1}^{\infty} a_n e^{-n^2\pi^2\nu t} n \sin(n\pi x)}{a_0 + \sum_{n=1}^{\infty} a_n e^{-n^2\pi^2\nu t} \cos(n\pi x)}
  $$
  係数（数値積分で算出）:
  $$
  a_0 = \int_0^1 \exp\left(-\frac{1-\cos(\pi x)}{2\pi\nu}\right) dx, \quad
  a_n = 2\int_0^1 \exp\left(-\frac{1-\cos(\pi x)}{2\pi\nu}\right)\cos(n\pi x) dx
  $$

### Problem 3: Cole-Hopf / Fourier Series (Different IC)

- **領域**: $0 \le x \le 1$
- **初期条件**:
  $$
  u(x,0) = 4x(1-x)
  $$
- **境界条件**: $u(0,t)=u(1,t)=0$
- **厳密解**: Problem 2 と同形式だが係数積分が以下となる:
  $$
  a_0 = \int_0^1 \exp\left(-\frac{x^2(3-2x)}{3\nu}\right) dx, \quad
  a_n = 2\int_0^1 \exp\left(-\frac{x^2(3-2x)}{3\nu}\right)\cos(n\pi x) dx
  $$

### Problem 4: Shock-like Solution ($t_0=1$)

- **領域**: $0 \le x \le 8$
- **初期条件**: $t_0=1$ における厳密解の値を与える
- **境界条件**: $u(0,t)=u(8,t)=0$
- **厳密解**:
  $$
  u(x,t) = \frac{x/t}{1 + \sqrt{t/t_0} \exp\left(\frac{x^2}{4\nu t}\right)}, \quad t_0 = \exp\left(\frac{1}{8\nu}\right)
  $$

### Problem 5: Piecewise IC / Inhomogeneous BC

- **領域**: $0 \le x \le 12$
- **初期条件**:
  $$
  u(x,0) = \begin{cases} 1 & (0 \le x \le 5) \\ 6-x & (5 < x \le 6) \\ 0 & (6 < x \le 12) \end{cases}
  $$
- **境界条件**:
  - $u(0,t) = 1$
  - $u(12,t) = 0$
- **厳密解**: なし（物理挙動の安定性を検証）

