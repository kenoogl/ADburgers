using Test
using ADburgers

# Define dummy boundary conditions for testing
dummy_bc_time(t) = 0.0
dummy_bc_coef(k, t) = 0.0

@testset "Taylor Engine Unit Tests" begin
    # Parameters
    N = 10
    h = 0.1
    x = collect(0:N) * h
    ν = 0.01
    K = 1 # Test order 1 generation (computes u_1 from u_0)
    
    # Consistent BCs for u=x^2 on [0,1]
    # left: u(0)=0. right: u(1)=1.
    # Taylor coeffs for BC:
    # left: (u(0))_0 = 0, others 0
    # right: (u(1))_0 = 1, others 0
    bc_left_coef(k, t) = 0.0
    bc_right_coef(k, t) = (k == 0 ? 1.0 : 0.0)
    
    # ProblemSpec mock
    prob = ProblemSpec(
        0.0, 1.0, ν, 0.0, 1.0,
        x -> x^2, 
        dummy_bc_time, dummy_bc_time,
        bc_left_coef, bc_right_coef
    )
    
    # Initialize arrays
    # U size is (N+1, K+1) => (11, 2)
    # Recall design: U[i+1, k+1]
    U = zeros(N+1, K+1)
    
    # Initialize u_0 = x^2 (internal points)
    # Boundary points will be overwritten by bc_coef in taylor_coeff!
    # But for this test, let's assume bc_coef consistent with 0 for simplicity or check internal only.
    # Actually, u_0 = x^2 implies u(0)=0 (left BC=0), u(1)=1 (right BC=1).
    # Let's set BC accordingly if we want consistency, but taylor_coeff! overwrites boundaries.
    # So we should check INTERNAL points.
    
    for i in 1:N+1
        U[i, 1] = x[i]^2
    end
    
    T1 = zeros(N+1)
    T2 = zeros(N+1)
    T3 = zeros(N+1)
    arrays = ADburgers.TaylorArrays(U, T1, T2, T3)
    
    # Run taylor_coeff!
    # We expect u_1 to be computed.
    t = 0.0
    ADburgers.taylor_coeff!(arrays, t, K, N, h, ν, prob)
    
    # Verification
    # Expected u_1 = -u u_x + ν u_xx
    # u = x^2
    # u_x = 2x
    # u_xx = 2
    # u_1 = -(x^2)(2x) + ν(2) = -2x^3 + 2ν
    
    # Check internal points i=2:N (indices 2:10)
    for iJ in 2:N
        xi = x[iJ]
        expected = -2 * xi^3 + 2 * ν
        computed = U[iJ, 2] # k=1 is index 2
        
        # Numerical differentiation error will exist.
        # Central difference for u_x of x^2 is EXACT for quadratic mesh?
        # ( (x+h)^2 - (x-h)^2 ) / 2h = (x^2 + 2xh + h^2 - (x^2 - 2xh + h^2)) / 2h = 4xh / 2h = 2x. EXACT.
        # Central difference for u_xx of x^2?
        # ( (x+h)^2 - 2x^2 + (x-h)^2 ) / h^2 = (2x^2 + 2h^2 - 2x^2) / h^2 = 2. EXACT.
        # So it should be close to floating point precision.
        
        @test isapprox(computed, expected, atol=1e-12)
    end
end
