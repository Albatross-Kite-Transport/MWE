using Interpolations, BenchmarkTools, StaticArrays

function test_interps()
    # Create test data
    n = 5  # Size of test matrix
    alphas = range(-10, 10, n)  # Angle of attack range: -10° to 10°
    d_trailing_edge_angles = range(0, 30, n)  # Trailing edge deflection: 0° to 30°

    # Create cl_matrix with a simple aerodynamic model
    cl_matrix = zeros(n, n)
    for (i, alpha) in enumerate(alphas)
        for (j, delta) in enumerate(d_trailing_edge_angles)
            # Simple linear combination of angle of attack and trailing edge deflection
            cl_matrix[i,j] = 2π * (deg2rad(alpha) + 0.3 * deg2rad(delta))
        end
    end

    cl_interp = extrapolate(scale(interpolate(cl_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)
    @time cl_interp(0.0, 0.0)
    @time cl_interp(0.0, 0.0)

    cl_interp2 = (α) -> cl_interp(α, 0.0) + cl_interp(α, 0.0)
    @time cl_interp2(0.0)
    @time cl_interp2(0.0)
end

test_interps()
