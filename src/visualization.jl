using RecipesBase

@recipe function f(sol::Solution; t_indices=nothing)
    xguide --> "x"
    yguide --> "u(x,t)"
    legend --> :best
    
    # Determine which time indices to plot
    # Default: Start, End, and a few middle ones (max 5 lines)
    if t_indices === nothing
        Nt = length(sol.t)
        if Nt <= 5
            inds = 1:Nt
        else
            inds = unique(round.(Int, range(1, Nt, length=5)))
        end
    else
        inds = t_indices
    end
    
    for i in inds
        t_val = sol.t[i]
        @series begin
            label --> "t = $(round(t_val, digits=3))"
            sol.x, sol.u[:, i]
        end
    end
end
