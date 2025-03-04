using Distributed
using Xfoil

const procs = addprocs()
println(test_x)
@everywhere begin
    using Xfoil
    function solve_alpha(x, y)
        Xfoil.set_coordinates(x, y)
        Xfoil.solve_alpha(0.0; mach=0.1)
        return nothing
    end
end

try
    x = collect(0:0.1:1)
    y = x .^ 2
    @show x y
    @sync @distributed for j in 1:10
        solve_alpha(x, y)
    end
catch e
    println(e)
finally
    println("removing processes")
    rmprocs(procs)
end