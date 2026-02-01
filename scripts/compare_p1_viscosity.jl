using ADburgers
using Plots
using Printf

function compare_viscosities()
    viscosities = [0.1, 0.2, 0.5, 1.0]
    N = 40
    t_max = 0.001
    dt = t_max
    K_range = 2:4
    
    # Create output directory
    mkpath("plots")

    for K in K_range
        p = plot(
            title="Problem 1: Viscosity Comparison (t=$(t_max), K=$(K))",
            xlabel="x",
            ylabel="u"
        )

        for ν in viscosities
            println("Running ν=$ν, K=$K...")
            # Use default a=2.0 as per Problem 1 definition
            prob_raw = get_problem(1; ν=ν, a=2.0)
            # Override tmax to match Figure 1 conditions
            prob = ProblemSpec(
                prob_raw.a, prob_raw.b, prob_raw.ν, prob_raw.t0, t_max,
                prob_raw.u0, prob_raw.bc_left, prob_raw.bc_right,
                prob_raw.bc_left_coeff, prob_raw.bc_right_coeff,
                prob_raw.exact_solution
            )
            
            save_times = [prob.t0, t_max]
            spec = SolverSpec(N, dt, K, save_times)
            
            try
                sol = solve(prob, spec)
                
                # Extract final solution
                u_final = sol.u[:, end]
                x = range(prob.a, prob.b, length=length(u_final))
                
                plot!(p, x, u_final, label="ν=$(ν)", lw=2)
                println("  Max u: $(maximum(u_final))")
            catch e
                println("  Failed for ν=$ν: $e")
            end
        end
        
        outfile = "plots/comparison_p1_viscosity_K$(K).png"
        savefig(p, outfile)
        println("Saved plot to $outfile")
    end
end

compare_viscosities()
