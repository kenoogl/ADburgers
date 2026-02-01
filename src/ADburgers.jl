module ADburgers

using Printf
using QuadGK

# Submodules / Files
include("types.jl")
include("core.jl")
include("problems.jl")
include("factory.jl")
include("analysis.jl")
include("visualization.jl")

# Exports
export ProblemSpec, SolverSpec, Solution, ReferenceSpec
export get_problem, solve

end # module
