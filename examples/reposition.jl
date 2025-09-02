# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using SymbolicAWEModels, KiteUtils, LinearAlgebra, ControlPlots, OrdinaryDiffEqBDF

# Initialize the model
set = Settings("system.yaml")
set.sample_freq = 100
sam = SymbolicAWEModel(set)
init!(sam)
find_steady_state!(sam)

# Run the simulation with repositioning
lg = SymbolicAWEModels.sim_reposition!(
    sam,
    total_time=10.0,
    reposition_interval_s=0.5,
    target_elevation_deg=50.0,
    target_azimuth_deg=10.0
)

# Plot the results
plot(sam.sys_struct, lg)


