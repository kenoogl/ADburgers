using LinearAlgebra
using Statistics

"""
    errors(sol::Solution, prob::ProblemSpec; t::Float64) -> Dict{String, Float64}

Computes L2 and L_inf errors against the exact solution at time t.
Returns a dictionary with keys "L2" and "Linf".
If exact_solution is missing, returns empty dict.
"""
function errors(sol::Solution, prob::ProblemSpec; t::Float64)
    if prob.exact_solution === nothing
        return Dict{String, Float64}()
    end

    # Find closest time index
    idx = findmin(abs.(sol.t .- t))[2]
    current_time = sol.t[idx]
    
    # Check if time matches closely enough
    if abs(current_time - t) > 1e-12
        @warn "Requested time $t not found in solution. Using closest time $current_time"
    end
    
    u_num = sol.u[:, idx]
    x = sol.x
    
    u_exact = [prob.exact_solution(xi, current_time) for xi in x]
    
    diff = u_num .- u_exact
    
    # L_inf = max(|diff|)
    l_inf = norm(diff, Inf)
    
    # L2 = sqrt( sum(diff^2) * h )  (Composite Trapezoidal or similar approximations)
    dx = x[2] - x[1]
    # Trapezoidal rule for integration:
    # int f^2 dx approx dx * (sum(f^2) - 0.5*f[1]^2 - 0.5*f[end]^2)
    sum_sq = sum(diff.^2) - 0.5 * diff[1]^2 - 0.5 * diff[end]^2
    l2 = sqrt(sum_sq * dx)
    
    return Dict("L2" => l2, "Linf" => l_inf)
end

"""
    check_acceptance(sol::Solution, prob_id::Int, prob::ProblemSpec) -> Dict{String, Bool}

Verifies if the solution meets the Acceptance Criteria (ACC) defined in requirements.
"""
function check_acceptance(sol::Solution, prob_id::Int, prob::ProblemSpec)
    results = Dict{String, Bool}()
    
    if prob_id == 1
        # ACC-P1-T1: |u_num - u_exact| <= 5e-7
        # Specified points? Usually max error over domain.
        res = errors(sol, prob, t=prob.tmax) # Check at final time usually
        if haskey(res, "Linf")
            results["ACC-P1-T1"] = res["Linf"] <= 5e-7
        else
            results["ACC-P1-T1"] = false
        end
        
    elseif prob_id == 4
        # ACC-P4-T8: Linf <= 1e-6
        t_eval = sol.t[end]
        res = errors(sol, prob, t=t_eval)
        if haskey(res, "Linf")
            results["ACC-P4-T8"] = res["Linf"] <= 1e-6
        else
            results["ACC-P4-T8"] = false
        end

    elseif prob_id == 5
        # ACC-P5-FIG: BC preservation, max overshoot, No NaN/Inf
        # 1. BC Preservation
        # Left BC u(0)=1, Right BC u(12)=0
        tol = 1e-12
        bc_ok = all(abs.(sol.u[1, :] .- 1.0) .< tol) && all(abs.(sol.u[end, :] .- 0.0) .< tol)
        
        # 2. Overshoot
        # max(u) <= 1 + 1e-3
        max_u = maximum(sol.u)
        overshoot_ok = max_u <= 1.0 + 1e-3
        
        # 3. No NaN/Inf (checked during solve usually, but double check)
        finite_ok = all(isfinite.(sol.u))
        
        results["ACC-P5-FIG"] = bc_ok && overshoot_ok && finite_ok
    end
    
    return results
end

"""
    ref_self_convergence(prob_id::Int, t::Float64, xs::Vector{Float64}, refspec::ReferenceSpec) -> Bool

Checks if the Fourier reference solution has converged by comparing N_fourier and 2*N_fourier.
Requirement: max difference <= 1e-10.
"""
function ref_self_convergence(prob_id::Int, t::Float64, xs::Vector{Float64}, refspec::ReferenceSpec)
    if prob_id == 2
        p_nf = ADburgers.problem2(n_max=refspec.Nfourier)
        p_2nf = ADburgers.problem2(n_max=2*refspec.Nfourier)
    elseif prob_id == 3
        p_nf = ADburgers.problem3(n_max=refspec.Nfourier)
        p_2nf = ADburgers.problem3(n_max=2*refspec.Nfourier)
    else
        return true # Not a Fourier problem
    end
    
    u_nf = [p_nf.exact_solution(x, t) for x in xs]
    u_2nf = [p_2nf.exact_solution(x, t) for x in xs]
    
    diff = norm(u_nf .- u_2nf, Inf)
    
    passed = diff <= refspec.self_conv_tol
    
    if !passed
        @warn "Reference solution self-convergence failed at t=$t. Diff=$diff > $(refspec.self_conv_tol)"
    end
    
    return passed
end
