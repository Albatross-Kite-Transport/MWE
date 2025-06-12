using CairoMakie

# Parameters
λ_d = 0.001  # Dangerous failure rate (failures per year)
T = 1.0      # Test interval (years)
n_contactors = 2  # Number of contactors for series/parallel examples
time = 0:0.01:2*T  # Time range over two test intervals

# PFD for a single contactor over time (reset at test interval T)
function pfd_single(t, λ_d, T)
    t_mod = mod(t, T)  # Reset at each test interval
    return λ_d * t_mod  # PFD increases linearly until reset
end

# PFD_avg calculation for a single contactor
pfd_avg_single = λ_d * T / 2

# PFD_avg for series arrangement (sum of individual PFDs)
pfd_avg_series = n_contactors * pfd_avg_single

# PFD_avg for parallel arrangement (product of individual PFDs)
pfd_avg_parallel = pfd_avg_single^n_contactors

# Data for PFD over time
pfd_values = [pfd_single(t, λ_d, T) for t in time]

# Plot 1: PFD over time for a single contactor
fig1 = Figure(resolution=(800, 400))
ax1 = Axis(fig1[1, 1],
    title="PFD over Time for a Single Contactor",
    xlabel="Time (years)",
    ylabel="PFD",
    xticks=0:0.5:2,
    yticks=0:0.0005:0.002)

lines!(ax1, time, pfd_values, color=:blue, label="PFD(t) = λ_d * t (reset at T)")
hlines!(ax1, [pfd_avg_single], color=:red, linestyle=:dash, label="PFD_avg = λ_d * T / 2")
axislegend(ax1, position=:lt)
ylims!(ax1, 0, 0.002)

# Save as SVG (vector graphics)
save("pfd_over_time.svg", fig1)

# Plot 2: Comparison of PFD_avg for series and parallel arrangements
fig2 = Figure(resolution=(800, 400))
ax2 = Axis(fig2[1, 1],
    title="PFD_avg for Series vs Parallel Arrangements",
    xlabel="Configuration",
    ylabel="PFD_avg",
    xticks=(1:3, ["Single", "Series (n=2)", "Parallel (n=2)"]))

# Bar plot data
configs = [pfd_avg_single, pfd_avg_series, pfd_avg_parallel]
barplot!(ax2, 1:3, configs, color=[:blue, :orange, :green])
text!(ax2, 1:3, configs .+ 0.00005, text=["0.0005", "0.001", "0.00000025"], align=(:center, :bottom))

# Adjust y-axis for visibility
ylims!(ax2, 0, maximum(configs) * 1.2)

# Save as SVG (vector graphics)
save("pfd_series_parallel.svg", fig2)

# Optionally display the plots (note: display may still rasterize in some environments)
display(fig1)
display(fig2)
