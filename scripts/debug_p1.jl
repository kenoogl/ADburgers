using ADburgers
using Printf
using LinearAlgebra

function test_p1_convergence()
    println("Debugging Problem 1 Convergence...")
    
    # N values to test
    Ns = [40, 80, 160]
    errors_list = []
    
    for N in Ns
        println("\nRunning N=$N...")
        # Manually define P1 with a=1.1
        # Need to copy-paste or modify source?
        # problems.jl hardcodes a=2.
        # Let's verify what happens if we could change 'a'.
        # Since I cannot pass 'a' to get_problem(1), I must edit problems.jl temporarily or define here.
        # Easier to edit problems.jl to accept keyword 'a'.
        prob = get_problem(1; a=1.01)
        dt = 0.0001 * (40/N) # Scale dt with dx for stability/accuracy parity
        steps = round(Int, 1.0/dt)
        save_times = [0.0, 1.0]
        
        solver_spec = SolverSpec(N, dt, 6, save_times)
        sol = solve(prob, solver_spec)
        
        errs = ADburgers.errors(sol, prob, t=1.0)
        push!(errors_list, errs["Linf"])
        
        max_u = maximum(sol.u)
        println("  Max u: $max_u")
        println("  Linf Error: $(errs["Linf"])")
    end
    
    println("\nConvergence Analysis:")
    for i in 1:(length(Ns)-1)
        ratio = errors_list[i] / errors_list[i+1]
        order = log2(ratio)
        println("N $(Ns[i]) -> $(Ns[i+1]): Ratio = $(round(ratio, digits=2)), Order = $(round(order, digits=2))")
    end
end

test_p1_convergence()
