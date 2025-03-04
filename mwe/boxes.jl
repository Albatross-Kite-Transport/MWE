using Interpolations
using Plots
using NonlinearSolve

# Read vertices from the .obj file
function read_vertices(filename)
    vertices = []
    open(filename) do file
        for line in eachline(file)
            if startswith(line, "v ") && !startswith(line, "vt") && !startswith(line, "vn")
                parts = split(line)
                x = parse(Float64, parts[2])
                y = parse(Float64, parts[3])
                z = parse(Float64, parts[4])
                push!(vertices, (x, y, z))
            end
        end
    end
    return vertices
end

function find_box_center_and_radius(vertices)
    r = zeros(2)
    v_min = zeros(3)
    v_tip = zeros(3)
    v_min .= Inf

    # find the vertex with smallest x in the middle of the kite
    for v in vertices
        if abs(v[2]) ≤ 0.1
            if v[1] < v_min[1]
                v_min .= v
            end
        end
    end

    # Find vertex furthest in y, -z direction
    max_score = -Inf
    v_tip .= 0.0
    for v in vertices
        # Score each vertex based on y and -z components
        # Higher y and lower z gives higher score
        score = v[2] - v[3]  # y - z
        if score > max_score
            max_score = score
            v_tip .= v
        end
    end

    function r_diff!(du, u, p)
        z = u[1]
        r .= Inf
        r[1] = sqrt(v_min[2]^2 + (v_min[3] - z)^2)
        r[2] = sqrt(v_tip[2]^2 + (v_tip[3] - z)^2)
        du[1] = r[1] - r[2]
        return nothing
    end

    z_min = Inf
    for v in vertices
        if abs(v[2]) ≤ 0.1
            if v[3] < z_min
                z_min = v[3]
            end
        end
    end

    prob = NonlinearProblem(r_diff!, [z_min], nothing)
    result = solve(prob, NewtonRaphson(; autodiff=AutoFiniteDiff(; relstep = 1e-3, absstep = 1e-3)); abstol = 1e-2)
    r_diff!(zeros(1), result, nothing)
    z = result[1]

    gamma_tip = atan(-v_tip[2], (v_tip[3] - z))
    @show rad2deg(gamma_tip)

    return z, r[1], gamma_tip
end

# Create interpolations for max and min x coordinates
function create_interpolations(vertices, box_center_z, radius, gamma_tip)
    gamma_range = range(-abs(gamma_tip)+1e-6, abs(gamma_tip)-1e-6, 100)
    vz_centered = [v[3] - box_center_z for v in vertices]
    
    max_xs = zeros(length(gamma_range))
    min_xs = zeros(length(gamma_range))
    areas  = zeros(length(gamma_range))
    
    for (j, gamma) in enumerate(gamma_range)
        max_xs[j] = -Inf
        min_xs[j] = Inf
        for (i, v) in enumerate(vertices)
            # Rotate y coordinate to check box containment
            rotated_y = v[2] * cos(gamma) - vz_centered[i] * sin(gamma)
            if gamma ≤ 0.0 && -0.5 ≤ rotated_y ≤ 0.0
                max_xs[j] = max(max_xs[j], v[1])
                min_xs[j] = min(min_xs[j], v[1])
            elseif gamma > 0.0 && 0.0 ≤ rotated_y ≤ 0.5
                max_xs[j] = max(max_xs[j], v[1])
                min_xs[j] = min(min_xs[j], v[1])
            end
        end
        area = abs(max_xs[j] - min_xs[j]) * gamma_range.step * radius
        last_area = j > 1 ? areas[j-1] : 0.0
        areas[j] = last_area + area
    end

    itp_max = linear_interpolation(gamma_range, max_xs)
    itp_min = linear_interpolation(gamma_range, min_xs)
    itp_area = linear_interpolation(gamma_range, areas)
    
    return (itp_max, itp_min, itp_area, gamma_range, max_xs, min_xs)
end

# Plot the interpolations
function plot_interpolations(itp_max, itp_min, itp_area, gammas, max_xs, min_xs)
    # Generate finer gamma values for smooth interpolation plots
    gamma_fine = range(minimum(gammas), maximum(gammas), length=1000)
    
    # Plot 1: Max/Min x coordinates
    p1 = plot(gamma_fine, itp_max.(gamma_fine),
        label="Maximum x",
        xlabel="Gamma (radians)",
        ylabel="x coordinate",
        title="Max and Min x Coordinates",
        linewidth=2
    )
    plot!(p1, gamma_fine, itp_min.(gamma_fine),
        label="Minimum x",
        linewidth=2
    )
    
    scatter!(p1, gammas, max_xs,
        label="",
        markercolor=:blue,
        markersize=3
    )
    scatter!(p1, gammas, min_xs,
        label="",
        markercolor=:red,
        markersize=3
    )
    
    # Plot 2: Area
    p2 = plot(gamma_fine, itp_area.(gamma_fine),
        label="Area",
        xlabel="Gamma (radians)",
        ylabel="Area",
        title="Cumulative Area",
        linewidth=2
    )
    
    # Combine plots in a layout
    plot(p1, p2, layout=(2,1), size=(800,600))
end

# Main execution
@time begin
    vertices = read_vertices("data/HL5_body.obj")
    box_center_z, radius, gamma_tip = find_box_center_and_radius(vertices)
    itp_max, itp_min, itp_area, gammas, max_xs, min_xs = create_interpolations(vertices, box_center_z, radius, gamma_tip)
    plot_interpolations(itp_max, itp_min, itp_area, gammas, max_xs, min_xs)
end