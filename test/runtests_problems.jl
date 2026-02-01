using Test
using ADburgers

@testset "Problem Definitions" begin
    # Test P1
    p1 = get_problem(1)
    @test p1.ν == 0.01
    @test p1.u0(0.5) ≈ (2*0.01*π*sin(0.5π))/(2 + cos(0.5π)) atol=1e-10
    if p1.exact_solution !== nothing
        @test p1.exact_solution(0.5, 0.0) ≈ p1.u0(0.5) atol=1e-10
    end

    # Test P4
    p4 = get_problem(4)
    @test p4.t0 == 1.0
    @test p4.exact_solution !== nothing
    # Check consistency of IC with exact solution at t0
    @test p4.u0(4.0) ≈ p4.exact_solution(4.0, 1.0) atol=1e-10

    # Test P5
    p5 = get_problem(5)
    @test p5.bc_left(0.0) == 1.0
    @test p5.u0(2.0) == 1.0
    @test p5.u0(5.5) == 0.5
    @test p5.u0(7.0) == 0.0
    @test p5.exact_solution === nothing

    # Test P2 (Integration check)
    p2 = get_problem(2)
    @test p2.u0(0.5) ≈ 1.0
    # Checking integration at t=0 against IC might be expensive or circular, 
    # but let's check it runs without error.
    if p2.exact_solution !== nothing
        val = p2.exact_solution(0.5, 0.0) 
        # t=0, exact solution should match IC = sin(pi*x)
        @test val ≈ sin(0.5π) atol=1e-3 # Fourier convergence tolerance
    end
end
