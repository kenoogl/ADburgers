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

function generate_p5_times(; ν=1.0, dx=0.05, dt=0.001, K=2, times=[0.0, 1.0, 2.0, 3.0, 4.0])
    mkpath("plots")

    prob_raw = get_problem(5; ν=ν)
    prob = override_tmax(prob_raw, maximum(times))

    domain_len = prob.b - prob.a
    N = Int(round(domain_len / dx))
    h = domain_len / N
    if abs(h - dx) > 1e-12
        error("dx does not evenly divide the domain. computed h=$(h), requested dx=$(dx)")
    end

    save_times = sort(times)
    if save_times[1] != prob.t0
        pushfirst!(save_times, prob.t0)
        save_times = unique(sort(save_times))
    end

    spec = SolverSpec(N, dt, K, save_times)
    sol = solve(prob, spec)

    # Diagnostics
    left_dev = maximum(abs.(sol.u[1, :] .- 1.0))
    right_dev = maximum(abs.(sol.u[end, :]))
    max_u = maximum(sol.u)
    println("Diagnostics: max_u=$(max_u), left_dev=$(left_dev), right_dev=$(right_dev)")
    for t in times
        idx = findfirst(==(t), sol.t)
        if idx !== nothing
            println("t=$(t): u(0)=$(sol.u[1, idx]), u(12)=$(sol.u[end, idx])")
        end
    end

    p = plot(
        title="Problem 5 (ν=$(ν), dx=$(dx), dt=$(dt), K=$(K))",
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

    outfile = "plots/Problem5_times.png"
    savefig(p, outfile)
    println("Saved plot to $outfile")
end

generate_p5_times()
