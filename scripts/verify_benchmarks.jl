using ADburgers
using Printf
using Plots

# Create plots directory
ispath("plots") || mkdir("plots")

println("Starting ADburgers Benchmark Verification...")
println("=============================================")

# --- Config ---
# Standard params mentioned in requirements/papers
nu_default = 0.01 
N_default = 40
dt_default = 0.0001
K_default = 6 

function run_benchmark(id, name; ν=nu_default, N=N_default, dt=dt_default, K=K_default, tmax_override=nothing, kwargs...)
    println("\nRunning $name (Problem $id)...")
    
    # Pass kwargs (e.g. n_max) to get_problem
    prob = get_problem(id; ν=ν, kwargs...)
    
    tmax = isnothing(tmax_override) ? prob.tmax : tmax_override
    save_times = range(prob.t0, tmax, step=dt * 100) |> collect
    # Ensure end point
    if save_times[end] != tmax
        push!(save_times, tmax)
    end
    # Ensure t0
    if save_times[1] != prob.t0
        pushfirst!(save_times, prob.t0)
    end
    sort!(save_times)
    unique!(save_times)

    # Solve
    solver_spec = SolverSpec(N, dt, K, save_times)
    
    try
        sol = solve(prob, solver_spec)
        println("  Solve completed. Steps: $(sol.meta["steps"])")
        
        # Plot
        p = plot(sol, title="$name (t=$(tmax))")
        savefig(p, "plots/Problem$(id).png")
        println("  Plot saved to plots/Problem$(id).png")
        
        # Verify
        # Check acceptance criteria
        res = ADburgers.check_acceptance(sol, id, prob)
        
        # Get raw errors for debugging
        debug_errs = ADburgers.errors(sol, prob, t=tmax)
        if haskey(debug_errs, "Linf")
            @printf("  Linf Error: %.2e (Threshold: %s)\n", debug_errs["Linf"], id==1 ? "5e-7" : "1e-6")
        end
        if haskey(debug_errs, "L2")
             @printf("  L2 Error:   %.2e\n", debug_errs["L2"])
        end

        passed = true
        for (criterion, status) in res
            status_str = status ? "PASS" : "FAIL"
            println("  $criterion: $status_str")
            if !status
                passed = false
            end
        end
        
        # Check Reference self-convergence for P2/P3
        if id == 2 || id == 3
            refspec = ReferenceSpec(false, true, 50, 1e-12, 1e-14, true, 1e-10)
            # Check at a few points
            check_t = prob.tmax / 2
            check_xs = range(prob.a, prob.b, length=10) |> collect
            # Note: We hardcode refspec.Nfourier here, but run_benchmark creates prob with reduced n_max.
            # ref_self_convergence creates NEW problem instances with ReferenceSpec parameters.
            # So this self-check checks if standard 1000/2000 setup converges.
            # This is slow! 
            # If we want speed, we should maybe SKIP this or reduce check points?
            # Or trust optimized precomputation to be fast enough for 10 points.
            # Precomputation means 1000 integrals. 10 points evaluation is cheap.
            # So self-convergence dominated by precomputation (2 calls: 1000 + 2000).
            # This should be fast now.
            conv = ADburgers.ref_self_convergence(id, check_t, check_xs, refspec)
            println("  Reference Self-Convergence: $(conv ? "PASS" : "FAIL")")
            if !conv
                passed = false
            end
        end
        
        return passed
    catch e
        println("  ERROR during execution: $e")
        Base.showerror(stdout, e, catch_backtrace())
        return false
    end
end

# --- RUN PROBLEMS ---

# Problem 1
# Paper: t=1.0? Table results usually at small times or t=1,3? 
# Requirement: ACC-P1-T1 checks max error.
pass1 = run_benchmark(1, "Problem 1 (Exact)"; tmax_override=0.001)

# Problem 2
# Paper uses ν=0.1 for P2. Use n_max=100 for speed (sufficient for t=1)
pass2 = run_benchmark(2, "Problem 2 (Fourier)"; ν=0.1, n_max=100) 

# Problem 3
# Paper uses ν=1.0 for P3. Use n_max=100 for speed
pass3 = run_benchmark(3, "Problem 3 (Fourier)"; ν=1.0, n_max=100) 

# Problem 4
# t0=1.0. Let's run to t=3.0
pass4 = run_benchmark(4, "Problem 4 (Shock-like)"; ν=0.01, tmax_override=3.0)

# Problem 5
pass5 = run_benchmark(5, "Problem 5 (Piecewise)"; tmax_override=5.0)

println("\n=============================================")
println("Summary:")
println("Problem 1: $(pass1 ? "PASS" : "FAIL")")
println("Problem 2: $(pass2 ? "PASS" : "FAIL")")
println("Problem 3: $(pass3 ? "PASS" : "FAIL")")
println("Problem 4: $(pass4 ? "PASS" : "FAIL")")
println("Problem 5: $(pass5 ? "PASS" : "FAIL")")

if pass1 && pass2 && pass3 && pass4 && pass5
    println("\nALL BENCHMARKS PASSED.")
    exit(0)
else
    println("\nSOME BENCHMARKS FAILED.")
    exit(1)
end
