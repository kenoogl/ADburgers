using Test

@testset "ADburgers Tests" begin
    include("runtests_core.jl")
    include("runtests_problems.jl")
end
