using ADburgers
using Plots

function override_tmax(prob::ProblemSpec, tmax::Float64)
    return ProblemSpec(
        prob.a, prob.b, prob.ν, prob.t0, tmax,
        prob.u0, prob.bc_left, prob.bc_right,
        prob.bc_left_coeff, prob.bc_right_coeff,
        prob.exact_solution
    )
end

function generate_figure2()
    ν = 0.01
    N = 80
    dt = 0.0001
    K = 2
    times = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0]

    mkpath("plots")

    # Build Problem 2 without precomputing Cole-Hopf coefficients
    # (exact solution is not required for this figure)
    L = 0.0
    R = 1.0
    t0 = 0.0
    tmax = maximum(times)
    u0_func(x) = sin(π * x)
    bc_left(t) = 0.0
    bc_right(t) = 0.0
    bc_left_coef(k, t) = 0.0
    bc_right_coef(k, t) = 0.0
    prob = ProblemSpec(L, R, ν, t0, tmax, u0_func, bc_left, bc_right, bc_left_coef, bc_right_coef, nothing)

    save_times = sort(times)
    if save_times[1] != prob.t0
        pushfirst!(save_times, prob.t0)
        save_times = unique(sort(save_times))
    end

    spec = SolverSpec(N, dt, K, save_times)
    sol = solve(prob, spec)

    # Boundary diagnostics
    left_dev = maximum(abs.(sol.u[1, :]))
    right_dev = maximum(abs.(sol.u[end, :]))
    println("Max boundary deviation: left=$(left_dev), right=$(right_dev)")
    for t in times
        idx = findfirst(==(t), sol.t)
        if idx !== nothing
            println("t=$(t): u(0)=$(sol.u[1, idx]), u(1)=$(sol.u[end, idx])")
        end
    end

    p = plot(
        title="Problem 2: Figure 2 (ν=$(ν), N=$(N), K=$(K))",
        xlabel="x",
        ylabel="u(x,t)",
        legend=:best
    )

    for t in times
        idx = findfirst(==(t), sol.t)
        if idx === nothing
            error("Requested time $t not found in solution output.")
        end
        plot!(p, sol.x, sol.u[:, idx], label="t=$(t)")
    end

    outfile = "plots/Problem2_Figure2.png"
    savefig(p, outfile)
    println("Saved plot to $outfile")
end

generate_figure2()
