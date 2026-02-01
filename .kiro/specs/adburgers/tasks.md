# 実装タスク: adburgers

## フェーズ 1: コア実装 (v0.1.0)
- [x] **初期セットアップ** <!-- id: 100 -->
    - [x] パッケージ構造の作成 (`Pkg.generate`) <!-- id: 101 -->
    - [x] 型定義の実装 `types.jl` (`ProblemSpec`, `SolverSpec`) <!-- id: 102 -->
    - [x] `Core.jl` スケルトンの実装 (空関数) <!-- id: 103 -->

- [x] **Taylor エンジン (コアロジック)** <!-- id: 200 -->
    - [x] `taylor_coeff!` の実装 (再帰式) <!-- id: 201 -->
        - [x] 単体テスト検証: 単純多項式 `t*x^2`
    - [x] `horner_update!` の実装 (時間進行) <!-- id: 202 -->
    - [x] `solve` の実装 (ループ, BC, 保存処理) <!-- id: 203 -->

- [x] **問題定義 & Factory** <!-- id: 300 -->
    - [x] `problems.jl`: Problem 1 実装 (厳密解) <!-- id: 301 -->
    - [x] `problems.jl`: Problem 4 実装 (区分的/厳密解) <!-- id: 302 -->
    - [x] `problems.jl`: Problem 5 実装 (衝撃波) <!-- id: 303 -->
    - [x] `factory.jl`: `get_problem(id)` の実装 <!-- id: 304 -->

## フェーズ 2: 検証・解析 (v0.2.0)
- [x] **参照解 (Problem 2/3)** <!-- id: 400 -->
    - [x] `ReferenceSpec` 構造体の実装 <!-- id: 401 -->
    - [x] `Cole-Hopf` Fourier 積分器の実装 (`QuadGK`) <!-- id: 402 -->
    - [x] `ref_self_convergence` (自己収束判定) の実装 <!-- id: 403 -->
    - [x] `problems.jl`: Problem 2 & 3 実装 <!-- id: 404 -->

- [/] **解析 & 可視化** <!-- id: 500 -->
    - [x] `analysis.jl` の実装 (L-inf, L2 誤差) <!-- id: 501 -->
    - [x] `visualization.jl` の実装 (RecipesBase/Plots) <!-- id: 502 -->
    - [x] `scripts/verify_benchmarks.jl` の作成 <!-- id: 503 -->

## フェーズ 3: 受入テスト・ドキュメント (v1.0.0)
- [/] **最終検証** <!-- id: 600 -->
    - [x] ACC-P1-T1 検証 (P1 Exact) <!-- id: 601 -->
    - [x] ACC-P4-T8 検証 (P4 Exact) <!-- id: 602 -->
    - [/] ACC-P5-FIG 検証 (P5 Qualitative) <!-- id: 603 -->
    - [x] ACC-P2/P3 検証 (Fourier Ref) <!-- id: 604 -->

- [ ] **ドキュメント** <!-- id: 700 -->
    - [ ] `README.md` の完成 <!-- id: 701 -->
    - [ ] API ドキュメント生成 (任意) <!-- id: 702 -->
