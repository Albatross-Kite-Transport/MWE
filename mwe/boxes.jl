using Interpolations
using Plots

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

# Create interpolations for max and min x coordinates
function create_interpolations(vertices; gamma_range=-deg2rad(83.7):0.01:deg2rad(83.7), box_center_z=-2.45)
    vz_centered = [v[3] - box_center_z for v in vertices]
    
    max_xs = Float64[]
    min_xs = Float64[]
    areas  = Float64[]
    gammas = Float64[]
    
    for gamma in gamma_range
        cos_g = cos(gamma)
        sin_g = sin(gamma)
        x_inside = Float64[]
        
        for (i, v) in enumerate(vertices)
            # Rotate y coordinate to check box containment
            rotated_y = v[2] * cos_g + vz_centered[i] * sin_g
            if abs(rotated_y) ≤ 0.1
                push!(x_inside, v[1])
            end
        end
        
        if !isempty(x_inside)
            push!(gammas, gamma)
            max_x = maximum(x_inside)
            min_x = minimum(x_inside)
            push!(max_xs, max_x)
            push!(min_xs, min_x)
            area = abs(max_x - min_x) * 0.01 * 1.63
            last_area = length(areas) > 0 ? areas[end] : 0.
            push!(areas, area + last_area)
        end
    end
    
    itp_max = linear_interpolation(gammas, max_xs)
    itp_min = linear_interpolation(gammas, min_xs)
    itp_area = linear_interpolation(gammas, areas)
    
    return (itp_max, itp_min, itp_area, gammas, max_xs, min_xs)
end

# Plot the interpolations
function plot_interpolations(itp_max, itp_min, itp_area, gammas, max_xs, min_xs)
    # Generate finer gamma values for smooth interpolation plots
    gamma_fine = range(minimum(gammas), maximum(gammas), length=1000)
    
    # Evaluate interpolations at finer gamma values
    max_x_fine = itp_max.(gamma_fine)
    min_x_fine = itp_min.(gamma_fine)
    
    # Create the plot
    plot(gamma_fine, max_x_fine, label="Maximum x", xlabel="Gamma (radians)", ylabel="x coordinate", title="Interpolated Max and Min x Coordinates", linewidth=2)
    plot!(gamma_fine, min_x_fine, label="Minimum x", linewidth=2)
    
    # Overlay the original data points
    scatter!(gammas, max_xs, label="", markercolor=:blue, markersize=3)
    scatter!(gammas, min_xs, label="", markercolor=:red, markersize=3)
    
    area_fine = itp_area.(gamma_fine)
    plot!(gamma_fine, area_fine, label="Area", xlabel="Gamma (radians)", ylabel="x coordinate", title="Interpolated Area", linewidth=2)
    display(plot!())
end

# Main execution
@time begin
    vertices = read_vertices("kite.obj")
    itp_max, itp_min, itp_area, gammas, max_xs, min_xs = create_interpolations(vertices)
    plot_interpolations(itp_max, itp_min, itp_area, gammas, max_xs, min_xs)
end