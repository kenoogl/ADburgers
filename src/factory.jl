# factory.jl

"""
    get_problem(id::Int; kwargs...) -> ProblemSpec

Returns the problem specification for the given ID (1-5).
"""
function get_problem(id::Int; kwargs...)
    if id == 1
        return problem1(; kwargs...)
    elseif id == 2
        return problem2(; kwargs...)
    elseif id == 3
        return problem3(; kwargs...)
    elseif id == 4
        return problem4(; kwargs...)
    elseif id == 5
        return problem5(; kwargs...)
    else
        throw(ArgumentError("Problem ID must be between 1 and 5"))
    end
end

