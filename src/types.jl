# types.jl

"""
    BCTime = Function

Signature: `(t::Float64) -> Float64`
Describes time-dependent boundary conditions.
"""
const BCTime = Function

"""
    BCCoef = Function

Signature: `(k::Int, t::Float64) -> Float64`
Describes Taylor coefficients of boundary conditions.
"""
const BCCoef = Function

"""
    ProblemSpec

Defines the physical problem parameters and conditions.
"""
struct ProblemSpec
    a::Float64
    b::Float64
    ν::Float64
    t0::Float64
    tmax::Float64
    
    u0::Function            # (x::Float64) -> Float64
    bc_left::BCTime
    bc_right::BCTime
    bc_left_coeff::BCCoef
    bc_right_coeff::BCCoef
end

"""
    SolverSpec

Defines the numerical solver parameters.
"""
struct SolverSpec
    N::Int
    Δt::Float64
    K::Int
    save_times::Vector{Float64}  # Must align with time grid
end

"""
    Algorithm: Taylor Arrays

Pre-allocated arrays for Taylor series recursion.
"""
struct TaylorArrays
    U::Matrix{Float64}  # (N+1, K+1)
    T1::Vector{Float64} # (N+1)
    T2::Vector{Float64} # (N+1)
    T3::Vector{Float64} # (N+1)
end

"""
    ReferenceSpec

Parameters for generating and verifying against reference solutions (Cole-Hopf/Fourier).
"""
struct ReferenceSpec
    has_exact::Bool
    has_reference::Bool
    Nfourier::Int             # default 1000, min 500
    rtol::Float64             # default 1e-12
    atol::Float64             # default 1e-14
    enforce_self_convergence::Bool
    self_conv_tol::Float64    # 1e-10
end

# Default constructor for ReferenceSpec
function ReferenceSpec(;
    has_exact=false,
    has_reference=false,
    Nfourier=1000,
    rtol=1e-12,
    atol=1e-14,
    enforce_self_convergence=true,
    self_conv_tol=1e-10
)
    ReferenceSpec(has_exact, has_reference, Nfourier, rtol, atol, enforce_self_convergence, self_conv_tol)
end

"""
    Solution

Results of the simulation.
"""
struct Solution
    x::Vector{Float64}
    t::Vector{Float64}
    u::Matrix{Float64}         # size (N+1, Nt)
    meta::Dict{String,Any}
end
