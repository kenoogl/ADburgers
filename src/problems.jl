# problems.jl

"""
    problem1(; ν=0.01, N=40) -> ProblemSpec

**Problem 1**: Exact solution available.
Domain: [0, 1]
BC: u(0,t) = u(1,t) = 0
IC: u(x,0) = 2νπ sin(πx) / (a + cos(πx)), a=2 (>1)
Exact: Eq (13) in paper.
"""
function problem1(; ν=0.01, N=40) # N argument is mainly for solver, but kept here if needed for defaults context
    a_param = 2.0 # Fixed as per problem description "a=2" or "a>1"
    
    L = 0.0
    R = 1.0
    t0 = 0.0
    tmax = 1.0 # Default, can be overridden by SolverSpec usually, but ProblemSpec defines physics
    
    # IC
    u0_func(x) = (2 * ν * π * sin(π * x)) / (a_param + cos(π * x))
    
    # BCs (Homogeneous Dirichlet)
    bc_left(t) = 0.0
    bc_right(t) = 0.0
    
    # Taylor coeffs for BC (all 0 for homogeneous constant)
    bc_left_coef(k, t) = 0.0
    bc_right_coef(k, t) = 0.0
    
    # Exact Solution
    function exact_sol(x, t)
        exp_term = exp(-π^2 * ν * t)
        numer = 2 * ν * π * exp_term * sin(π * x)
        denom = a_param + exp_term * cos(π * x)
        return numer / denom
    end
    
    return ProblemSpec(L, R, ν, t0, tmax, u0_func, bc_left, bc_right, bc_left_coef, bc_right_coef, exact_sol)
end


"""
    problem4(; ν=0.01) -> ProblemSpec

**Problem 4**: Shock-like solution.
Domain: [0, 8]
Time: Starts at t=1.0
BC: u(0,t) = u(8,t) = 0
IC (at t=1): Derived from exact solution Eq (16)
"""
function problem4(; ν=0.01)
    L = 0.0
    R = 8.0
    t0 = 1.0 # Starts at t=1
    tmax = 3.0 # Arbitrary default
    
    t0_param = exp(1.0 / (8.0 * ν))
    
    # Exact Solution Eq (16)
    function exact_sol(x, t)
        if t < 1.0
            error("Problem 4 defined for t >= 1")
        end
        # u(x,t) = (x/t) / (1 + sqrt(t/t0)*exp(x^2/(4νt)))
        term1 = x / t
        term2 = sqrt(t / t0_param)
        term3 = exp(x^2 / (4.0 * ν * t))
        return term1 / (1.0 + term2 * term3)
    end

    # IC at t0=1
    u0_func(x) = exact_sol(x, t0)
    
    # BCs (Homogeneous Dirichlet)
    bc_left(t) = 0.0
    bc_right(t) = 0.0
    bc_left_coef(k, t) = 0.0
    bc_right_coef(k, t) = 0.0
    
    return ProblemSpec(L, R, ν, t0, tmax, u0_func, bc_left, bc_right, bc_left_coef, bc_right_coef, exact_sol)
end

"""
    problem5(; ν=0.01) -> ProblemSpec

**Problem 5**: Piecewise IC, Non-homogeneous BC.
Domain: [0, 12]
BC: u(0,t) = 1, u(12,t) = 0
IC: Piecewise function
Exact: None (usually)
"""
function problem5(; ν=0.01)
    L = 0.0
    R = 12.0
    t0 = 0.0
    tmax = 5.0
    
    # Piecewise IC
    function u0_func(x)
        if 0.0 <= x <= 5.0
            return 1.0
        elseif 5.0 < x <= 6.0
            return 6.0 - x
        else # 6 < x <= 12
            return 0.0
        end
    end
    
    # BCs
    bc_left(t) = 1.0
    bc_right(t) = 0.0
    
    # Taylor coeffs for BC
    # Left is constant 1.0 => k=0 is 1.0, k>0 is 0.0
    # Right is constant 0.0 => all 0.0
    bc_left_coef(k, t) = (k == 0 ? 1.0 : 0.0)
    bc_right_coef(k, t) = 0.0
    

    return ProblemSpec(L, R, ν, t0, tmax, u0_func, bc_left, bc_right, bc_left_coef, bc_right_coef, nothing)
end

"""
    calc_cole_hopf_exact(x, t, ν, a0_integrand, an_integrand; n_max=1000)

Helper to calculate exact solution using Cole-Hopf Fourier Transformation.
"""
function calc_cole_hopf_exact(x, t, ν, a0_integrand, an_integrand; n_max=1000)
    # Compute a0
    # a0 = ∫_0^1 a0_integrand(x) dx
    a0, _ = quadgk(a0_integrand, 0.0, 1.0; rtol=1e-12)
    
    sum_numer = 0.0
    sum_denom = 0.0
    
    for n in 1:n_max
        # Compute an
        # an = 2 * ∫_0^1 an_integrand(x, n) dx
        # Optimization: Pass n to integrand closure if possible to avoid recompiling?
        # quadgk accepts function.
        an_func(y) = an_integrand(y, n)
        an_val, _ = quadgk(an_func, 0.0, 1.0; rtol=1e-12)
        an = 2.0 * an_val
        
        # Terms
        exp_factor = exp(-n^2 * π^2 * ν * t)
        
        sum_numer += an * exp_factor * n * sin(n * π * x)
        sum_denom += an * exp_factor * cos(n * π * x)
    end
    
    numer = 2 * π * ν * sum_numer
    denom = a0 + sum_denom
    
    return numer / denom
end

"""
    problem2(; ν=0.1, n_max=1000) -> ProblemSpec

**Problem 2**: Cole-Hopf solution.
IC: u(x,0) = sin(πx)
"""
function problem2(; ν=0.1, n_max=1000)
    L = 0.0
    R = 1.0
    t0 = 0.0
    tmax = 1.0
    
    u0_func(x) = sin(π * x)
    
    bc_left(t) = 0.0
    bc_right(t) = 0.0
    bc_left_coef(k, t) = 0.0
    bc_right_coef(k, t) = 0.0
    
    # Integrands for coefficients
    # Eq (15): exp(-(1-cos(πx))/(2πν))
    a0_int(x) = exp(-(1.0 - cos(π * x)) / (2.0 * π * ν))
    an_int(x, n) = exp(-(1.0 - cos(π * x)) / (2.0 * π * ν)) * cos(n * π * x)
    
    function exact_sol(x, t)
        return calc_cole_hopf_exact(x, t, ν, a0_int, an_int; n_max=n_max)
    end
    
    return ProblemSpec(L, R, ν, t0, tmax, u0_func, bc_left, bc_right, bc_left_coef, bc_right_coef, exact_sol)
end

"""
    problem3(; ν=1.0, n_max=1000) -> ProblemSpec

**Problem 3**: Cole-Hopf solution.
IC: u(x,0) = 4x(1-x)
"""
function problem3(; ν=1.0, n_max=1000)
    L = 0.0
    R = 1.0
    t0 = 0.0
    tmax = 1.0
    
    u0_func(x) = 4.0 * x * (1.0 - x)
    
    bc_left(t) = 0.0
    bc_right(t) = 0.0
    bc_left_coef(k, t) = 0.0
    bc_right_coef(k, t) = 0.0
    
    # Integrands for coefficients
    # exp(-x^2(3-2x)/(3ν))
    gamma(x) = (x^2 * (3.0 - 2.0 * x)) / (3.0 * ν)
    a0_int(x) = exp(-gamma(x))
    an_int(x, n) = exp(-gamma(x)) * cos(n * π * x)
    
    function exact_sol(x, t)
        return calc_cole_hopf_exact(x, t, ν, a0_int, an_int; n_max=n_max)
    end
    
    return ProblemSpec(L, R, ν, t0, tmax, u0_func, bc_left, bc_right, bc_left_coef, bc_right_coef, exact_sol)
end



