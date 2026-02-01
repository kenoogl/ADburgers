# core.jl

"""
    solve(prob::ProblemSpec, specs::SolverSpec) -> Solution

Main solver function.
Orchestrates the time stepping, data saving, and stability checks.
"""
function solve(prob::ProblemSpec, specs::SolverSpec)
    # Unpack specs
    N = specs.N
    dt = specs.Δt
    K = specs.K
    save_times = specs.save_times
    
    t0 = prob.t0
    tmax = prob.tmax
    h = (prob.b - prob.a) / N
    ν = prob.ν
    
    # Validation: save_times
    # Check if save_times align with dt grid
    for st in save_times
        m = (st - t0) / dt
        tol_m = 1000 * eps(Float64) * max(1.0, abs(m))
        if abs(m - round(m)) > tol_m
            throw(ArgumentError("save_time $st does not align with time grid (dt=$dt, t0=$t0)"))
        end
        if st < t0 || st > tmax
            throw(ArgumentError("save_time $st out of range [$t0, $tmax]"))
        end
    end
    # Ensure t0 is in save_times (design doc Requirement 3.2: must include t0)
    if !any(t -> isapprox(t, t0, atol=1e-12), save_times)
         throw(ArgumentError("save_times must include t0 ($t0)"))
    end

    # Initialize Memory
    U = zeros(N+1, K+1)
    T1 = zeros(N+1)
    T2 = zeros(N+1)
    T3 = zeros(N+1)
    arrays = TaylorArrays(U, T1, T2, T3)
    
    # Initialize Solution storage
    # Pre-allocate output u matrix
    Nt = length(save_times)
    u_out = zeros(N+1, Nt)
    
    # Initialize u0
    x_coords = range(prob.a, prob.b, length=N+1) |> collect
    for i in 1:N+1
        U[i, 1] = prob.u0(x_coords[i])
    end
    
    # Check BC consistency at t0
    # Tolerance 1e-12 as per 3.3
    tol_bc = 1e-12
    bc_l_t0 = prob.bc_left(t0)
    bc_r_t0 = prob.bc_right(t0)
    
    if abs(U[1, 1] - bc_l_t0) > tol_bc
        @warn "Initial condition at left boundary inconsistent with BC. Overwriting u0(a)."
        U[1, 1] = bc_l_t0
    end
    if abs(U[N+1, 1] - bc_r_t0) > tol_bc
        @warn "Initial condition at right boundary inconsistent with BC. Overwriting u0(b)."
        U[N+1, 1] = bc_r_t0
    end
    
    # Save initial state if t0 in save_times (checked above)
    # Find index of t0 in sorted save_times (assuming sorted? Design doc assumes sorted? Reco sorted)
    # Optimization: map time to index if sorted. We iterate.
    
    save_idx = 1
    # Assuming save_times is sorted. If strict, we sort it? SolverSpec inputs should be sorted ideally.
    if !issorted(save_times)
        error("save_times must be sorted")
    end
    
    if abs(save_times[1] - t0) < 1e-12
        u_out[:, 1] .= U[:, 1]
        save_idx += 1
    end
    
    # Time Loop
    t = t0
    step_count = 0
    max_steps = round(Int, (tmax - t0) / dt)
    
    cfl_max = 0.0
    diff_max = 0.0
    
    # Diffusion number constant: ν * dt / h^2
    diff_num = ν * dt / (h^2)
    
    while step_count < max_steps
        # Taylor Coefficients -> U[:, 2:K+1] computed
        taylor_coeff!(arrays, t, K, N, h, ν, prob)
        
        # Horner Update -> U[:, 1] updated to t+dt
        t_next = t + dt
        horner_update!(arrays, dt, K, N, prob.bc_left, prob.bc_right, t_next)
        
        t = t_next
        step_count += 1
        
        # Stability Check & NaN check
        # CFL = max|u| * dt / h
        # Diffusion = ν * dt / h^2 (Constant)
        
        current_max_u = maximum(abs, @view U[:, 1])
        cfl = current_max_u * dt / h
        
        cfl_max = max(cfl_max, cfl)
        diff_max = max(diff_max, diff_num)
        
        if !isfinite(current_max_u)
            error("NaN/Inf detected at t=$t")
        end
        
        # Save if needed
        # Check against save_times[save_idx]
        if save_idx <= Nt && abs(t - save_times[save_idx]) < 1e-12
            u_out[:, save_idx] .= U[:, 1]
            save_idx += 1
        end
    end
    
    meta = Dict{String, Any}(
        "CFL_max" => cfl_max,
        "diffusion_max" => diff_max,
        "steps" => step_count,
        "final_time" => t
    )
    
    return Solution(x_coords, save_times, u_out, meta)
end

"""
    step!(arrays::TaylorArrays, t::Float64, Δt::Float64, K::Int, N::Int, h::Float64, ν::Float64, prob::ProblemSpec)

Performs one time step update.
"""
function step!(arrays::TaylorArrays, t::Float64, Δt::Float64, K::Int, N::Int, h::Float64, ν::Float64, prob::ProblemSpec)
    # TODO: Implement step logic
    # 1. taylor_coeff!
    # 2. horner_update!
    error("Not implemented")
end

"""
    taylor_coeff!(arrays::TaylorArrays, t::Float64, K::Int, N::Int, h::Float64, ν::Float64, prob::ProblemSpec)

Computes Taylor coefficients (u_i)_k recursively for k=0 to K-1.
Result is stored in arrays.U.
"""
function taylor_coeff!(arrays::TaylorArrays, t::Float64, K::Int, N::Int, h::Float64, ν::Float64, prob::ProblemSpec)
    U = arrays.U
    T1 = arrays.T1
    T2 = arrays.T2
    T3 = arrays.T3
    
    # Pre-compute constants
    inv_2h = 1.0 / (2.0 * h)
    inv_h2 = 1.0 / (h^2)

    # Recursive calculation loop
    for k in 0:(K-1)
        # 1. Set boundary coefficients for order k
        # Julia indices: 1 (left), N+1 (right)
        U[1, k+1]     = prob.bc_left_coeff(k, t)
        U[N+1, k+1]   = prob.bc_right_coeff(k, t)
        
        # 2. Compute T1 (Advection part 1) and T3 (Diffusion) for current k
        # Loop over internal points iJ = 2:N
        @inbounds for iJ in 2:N
            u_next = U[iJ+1, k+1]
            u_prev = U[iJ-1, k+1]
            u_curr = U[iJ, k+1]
            
            T1[iJ] = (u_next - u_prev) * inv_2h
            T3[iJ] = (u_prev - 2.0 * u_curr + u_next) * inv_h2
        end
        
        # 3. Compute T2 (Convolution)
        # (T2_i)_k = sum_{j=0}^k (u_i)_j * (T1_i)_{k-j}
        # We need to re-compute T1_{k-j} on the fly or store it.
        # Design choice: Low allocation, re-compute T1 needed for convolution.
        
        # Initialize T2 with 0
        @views fill!(T2[2:N], 0.0)
        
        for j in 0:k
            # Calculate T1_{k-j}
            order_idx = (k - j) + 1 # 1-based index for U
            
            # Inner spatial loop for convolution
            @inbounds for iJ in 2:N
                # T1_{k-j} at i
                # Note: U is (N+1, K+1)
                t1_val = (U[iJ+1, order_idx] - U[iJ-1, order_idx]) * inv_2h
                
                # Add to T2
                T2[iJ] += U[iJ, j+1] * t1_val
            end
        end
        
        # 4. Update next coefficient (u_{k+1})
        # (u_i)_{k+1} = ( -T2 + ν*T3 ) / (k+1)
        denominator = 1.0 / (k + 1)
        @inbounds for iJ in 2:N
            U[iJ, (k+1)+1] = (-T2[iJ] + ν * T3[iJ]) * denominator
        end
    end
end


"""
    horner_update!(arrays::TaylorArrays, Δt::Float64, K::Int, N::Int, bc_left::BCTime, bc_right::BCTime, t_next::Float64)

Updates u (U[:, 1]) to next time step t+Δt using Horner's method with Taylor coefficients.
Also enforces boundary conditions at t+Δt.
"""
function horner_update!(arrays::TaylorArrays, Δt::Float64, K::Int, N::Int, bc_left::BCTime, bc_right::BCTime, t_next::Float64)
    U = arrays.U
    
    # Horner's method for internal points (i=2:N)
    # u(t+dt) = u_0 + u_1*dt + u_2*dt^2 + ... + u_K*dt^K
    #         = u_0 + dt*(u_1 + dt*(u_2 + ...))
    
    @inbounds for iJ in 2:N
        s = U[iJ, (K)+1] # Start with highest order coefficient
        for k in (K-1):-1:0
            s = s * Δt + U[iJ, k+1]
        end
        U[iJ, 1] = s
    end
    
    # Enforce boundary conditions at t_next
    U[1, 1]     = bc_left(t_next)
    U[N+1, 1]   = bc_right(t_next)
end

