using NonlinearSolve

function find_chord_lenght(w, t, m)
    total_area = w * (t + m) / 2
    target_area = total_area / 4
    function area(x)
        return m * x - 0.5 * x * ((m - t) * x / w)
    end
    function equations!(F, x, p)
        F[1] = area(x[1]) - target_area
        F[2] = area(x[2]) - 3 * target_area
    end

    x0 = [w / 4, 3 * w / 4]
    prob = NonlinearProblem(equations!, x0, nothing)

    result = solve(prob, NewtonRaphson())

    c1, c2 = result.u
    return c1, c2
end

w = 10.0
t = 2.0
m = 6.0

@time x1, x3 = find_chord_lenght(w, t, m)
@time x1, x3 = find_chord_lenght(w, t, m)
println("Vertical line positions: x1 = $x1, x3 = $x3")
