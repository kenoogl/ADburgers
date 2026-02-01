using ADburgers
using Printf
using TOML

const DEFAULT_SPEC = joinpath(@__DIR__, "p1_table_specs.toml")

function override_tmax(prob::ProblemSpec, tmax::Float64)
    return ProblemSpec(
        prob.a, prob.b, prob.ν, prob.t0, tmax,
        prob.u0, prob.bc_left, prob.bc_right,
        prob.bc_left_coeff, prob.bc_right_coeff,
        prob.exact_solution
    )
end

function sample_values(xs::Vector{Float64}, xgrid::Vector{Float64}, u::Vector{Float64}; tol=1e-12)
    vals = similar(xs)
    used_interp = false

    for (i, xi) in pairs(xs)
        # Clamp to domain endpoints
        if xi <= xgrid[1]
            vals[i] = u[1]
            continue
        elseif xi >= xgrid[end]
            vals[i] = u[end]
            continue
        end

        idx = searchsortedfirst(xgrid, xi)
        if idx <= length(xgrid) && abs(xgrid[idx] - xi) <= tol
            vals[i] = u[idx]
            continue
        elseif idx > 1 && abs(xgrid[idx - 1] - xi) <= tol
            vals[i] = u[idx - 1]
            continue
        end

        # Linear interpolation for off-grid points
        used_interp = true
        x0 = xgrid[idx - 1]
        x1 = xgrid[idx]
        u0 = u[idx - 1]
        u1 = u[idx]
        vals[i] = u0 + (u1 - u0) * (xi - x0) / (x1 - x0)
    end

    return vals, used_interp
end

function solve_problem1(ν::Float64, a::Float64, N::Int, dt::Float64, K::Int, t_eval::Float64)
    prob_raw = ADburgers.problem1(ν=ν, a=a)
    prob = override_tmax(prob_raw, t_eval)
    save_times = [prob.t0, t_eval]
    spec = SolverSpec(N, dt, K, save_times)
    sol = solve(prob, spec)

    u_num = sol.u[:, end]
    xgrid = sol.x
    u_exact = [prob.exact_solution(xi, t_eval) for xi in xgrid]

    return prob, xgrid, u_num, u_exact
end

function write_csv(path::String, xs::Vector{Float64}, u_num::Vector{Float64}, u_exact::Vector{Float64}, abs_err::Vector{Float64})
    open(path, "w") do io
        println(io, "x,u_num,u_exact,abs_err")
        for i in eachindex(xs)
            @printf(io, "%.16g,%.16g,%.16g,%.16g\n", xs[i], u_num[i], u_exact[i], abs_err[i])
        end
    end
end

function render_table(name::String, ν::Float64, a::Float64, N::Int, dt::Float64, K::Int, t_eval::Float64,
                      xs::Vector{Float64}, u_num::Vector{Float64}, u_exact::Vector{Float64}, abs_err::Vector{Float64};
                      digits::Int=6, x_digits::Int=4)
    println("== $(name) (ν=$(ν), a=$(a), N=$(N), dt=$(dt), K=$(K), t=$(t_eval)) ==")
    @printf("%10s %18s %18s %18s\n", "x", "u_num", "u_exact", "abs_err")
    for i in eachindex(xs)
        @printf("%10.*f %18.*f %18.*f %18.2e\n", x_digits, xs[i], digits, u_num[i], digits, u_exact[i], abs_err[i])
    end
    println()
end

function load_specs(path::String)
    if !isfile(path)
        error("Spec file not found: $path")
    end
    spec = TOML.parsefile(path)
    tables = get(spec, "tables", nothing)
    if tables === nothing || !(tables isa Vector)
        error("Spec must contain [[tables]] entries: $path")
    end
    return tables
end

function run()
    spec_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_SPEC
    tables = load_specs(spec_path)

    mkpath("tables")

    for tdef in tables
        name = String(get(tdef, "name", "Table"))
        ν_default = Float64(get(tdef, "nu", 0.0))
        nus = get(tdef, "nus", nothing)
        a = Float64(get(tdef, "a", 2.0))
        N = Int(tdef["N"])
        dt_default = Float64(get(tdef, "dt", 0.0))
        dts = get(tdef, "dts", nothing)
        t_eval = Float64(tdef["t"])
        xs = Float64.(tdef["x"])
        Ks = if haskey(tdef, "Ks")
            Int.(tdef["Ks"])
        elseif haskey(tdef, "K")
            [Int(tdef["K"])]
        else
            error("K or Ks must be provided for $name in $(spec_path).")
        end
        digits = Int(get(tdef, "digits", 6))
        x_digits = Int(get(tdef, "x_digits", 4))

        if isempty(xs)
            error("x list is empty for $name. Fill 'x' in $(spec_path).")
        end

        nu_list = nus === nothing ? [ν_default] : Float64.(nus)
        if isempty(nu_list) || nu_list[1] == 0.0
            error("nu or nus must be provided for $name in $(spec_path).")
        end

        dt_list = dts === nothing ? [dt_default] : Float64.(dts)
        if isempty(dt_list) || dt_list[1] == 0.0
            error("dt or dts must be provided for $name in $(spec_path).")
        end

        for ν in nu_list
            for dt in dt_list
                for K in Ks
                    prob, xgrid, u_num_grid, u_exact_grid = solve_problem1(ν, a, N, dt, K, t_eval)
                    u_num, used_interp = sample_values(xs, xgrid, u_num_grid)
                    u_exact, _ = sample_values(xs, xgrid, u_exact_grid)
                    abs_err = abs.(u_num .- u_exact)

                    if used_interp
                        @warn "Some x points are off-grid. Linear interpolation used." name K dt ν
                    end

                    render_table(name, ν, a, N, dt, K, t_eval, xs, u_num, u_exact, abs_err; digits=digits, x_digits=x_digits)

                    dt_tag = replace(@sprintf("%.8f", dt), "." => "p")
                    nu_tag = replace(@sprintf("%.4f", ν), "." => "p")
                    outfile = joinpath("tables", "$(name)_nu$(nu_tag)_K$(K)_dt$(dt_tag).csv")
                    write_csv(outfile, xs, u_num, u_exact, abs_err)
                    println("Saved: $outfile")
                end
            end
        end
    end
end

run()
