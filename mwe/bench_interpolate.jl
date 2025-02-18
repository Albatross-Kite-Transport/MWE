using Interpolations
using BenchmarkTools

# Test configurations
const DIMS = 1:6
const N_POINTS = 1_000
const GRID_POINTS = 5
const METHODS = [
    (BSpline(Cubic(Line(OnGrid()))), "Cubic"),
    (BSpline(Linear()), "Linear"),
    (BSpline(Quadratic(Line(OnGrid()))), "Quadratic")
]

# ...existing code for create_test_data...

function generate_test_points(dim::Int, grid_axes)
    # Generate evaluation points within grid bounds
    ε = 1e-10
    points = Vector{Vector{Float64}}(undef, N_POINTS)
    
    for i in 1:N_POINTS
        point = Vector{Float64}(undef, dim)
        for d in 1:dim
            # Generate point within bounds
            ax = grid_axes[d]
            point[d] = rand() * (last(ax) - first(ax) - 2ε) + first(ax) + ε
        end
        points[i] = point
    end
    
    return points
end

function evaluate_interpolation(itp, points)
    result = 0.0  # Just accumulate to avoid array allocation
    for point in points
        result += itp(point...)
    end
    return result
end

function benchmark_dimension(dim::Int)
    println("\nDimension $dim:")
    println("-" ^ 40)
    
    grid_axes, data = create_test_data(dim)
    test_points = generate_test_points(dim, grid_axes)
    
    for (method, name) in METHODS
        # Create interpolation with boundary conditions
        itp = interpolate(data, method)
        scaled_itp = scale(itp, grid_axes...)
        
        evaluate_interpolation(scaled_itp, test_points)
        
        # Benchmark just the interpolation evaluation
        b = @benchmark evaluate_interpolation($scaled_itp, $test_points)
        
        # Print results
        println("$name:")
        println("  Time (med):   $(median(b.times)/1e6) ms")
        println("  Memory:       $(b.memory/1024) KB")
        println("  Allocations:  $(b.allocs)")
    end
end

function run_benchmarks()
    println("Interpolation Benchmarks")
    println("======================")
    println("Configuration:")
    println("  Evaluation points: $N_POINTS")
    println("  Grid points:      $GRID_POINTS")
    
    for dim in DIMS
        benchmark_dimension(dim)
    end
end

# Run benchmarks
run_benchmarks()